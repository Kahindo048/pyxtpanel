"""How well does Stata agree with itself?  The level of agreement that closed forms can reach (B.3.4).

    python -m validation.self_agreement            # run Stata on a twin of the campaign, freeze
    python -m validation.self_agreement --check    # compare the frozen file with the golden data

The tolerances of the validation (``validation/compare.py``) say how far the package may be from
Stata.  Before they are tightened for the closed-form quantities (counts, moments, transformations)
one must know how far Stata is from *itself*, because an agreement better than that noise cannot
be asked of any implementation.  Three measures are taken on the 160 frozen scenarios of the
campaign (73 quantities each):

``stata_vs_twin``
    The same Stata battery on a *twin* of the data, in which the units of every scenario are
    renamed (the same labels, handed out in another order) and the rows shuffled.  ``xtset`` sorts
    by unit and period, so the unit order, and with it the order of the floating-point additions
    over units inside Stata, changes; every quantity is mathematically the same, except the ones
    that weigh the position of the unit (``*wsum``), which are left out of this comparison.  What
    remains is Stata's own numerical noise.  (Shuffling the rows alone proved to change nothing at
    all: ``xtset`` puts them back in the same order.)
``package_vs_stata``
    The package against the frozen Stata results (what ``tests/golden`` checks at 1e-9).
``display``
    What a person who reads the value on the screen instead of in ``r()`` has: the stored value
    written with ``%9.0g`` (nine significant digits), the width Stata's tables use, and ``%12.0g``,
    the default of ``string()``.

The error of a cell is ``|a - b| / max(|b|, 1)``: relative for the large values, absolute for the
small ones (a quantity that is zero in theory comes out as ``1e-17`` in one order and ``0`` in
another, which is no disagreement).  The result is written to ``tests/golden/stata_self_agreement
.json``; ``tests/golden/test_self_agreement.py`` checks it without Stata.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from validation import campaign
from validation.compare import read_matrix
from validation.registry import MATA_EXPORT

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "out"
GOLDEN = ROOT / "tests" / "golden" / "campaign_golden.json"
TARGET = ROOT / "tests" / "golden" / "stata_self_agreement.json"
TWIN_SEED = 20261119
DISPLAY_FORMATS = ("%9.0g", "%12.0g")


def error(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Per-cell error ``|a - b| / max(|b|, 1)``; ``inf`` where only one of the two is missing."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    out = np.abs(a - b) / np.maximum(np.abs(b), 1.0)
    out[np.isnan(a) & np.isnan(b)] = 0.0
    out[np.isnan(a) ^ np.isnan(b)] = np.inf
    return out


#: Quantities that weigh the position of the unit: they change when the units are renamed.
POSITION = [i for i, name in enumerate(campaign.COLUMNS) if name.endswith("wsum")]
UNIT_ORDER_FREE = [i for i in range(len(campaign.COLUMNS)) if i not in POSITION]


def twin() -> pd.DataFrame:
    """The frozen data with the units renamed in another order and the rows shuffled."""
    rng = np.random.default_rng(TWIN_SEED)
    stacked = pd.read_csv(campaign.DATA_FILE)
    parts = []
    for _, part in stacked.groupby("scenario", sort=True):
        labels = np.sort(part["id"].unique())
        renamed = dict(zip(labels, rng.permutation(labels), strict=True))
        part = part.assign(id=part["id"].map(renamed))
        parts.append(part.iloc[rng.permutation(len(part))])
    return pd.concat(parts, ignore_index=True)[["scenario", "id", "t", "x", "z"]]


def write_do() -> Path:
    """Write the twin data and the do-file that runs the battery on it and exports the cells."""
    OUT.mkdir(parents=True, exist_ok=True)
    data = twin()
    for column in ("scenario", "id", "t"):
        data[column] = data[column].astype("int32")
    dta = OUT / "sim_twin.dta"
    data.to_stata(dta, write_index=False, version=118)
    n_cols = len(campaign.COLUMNS)
    lines = [
        "version 17",
        "set more off",
        "clear all",
        MATA_EXPORT,
        campaign.MATA,
        campaign.stata_program(),
        f'use "{dta.as_posix()}", clear',
        f"mata: R = J(0, {n_cols}, .)",
        *[f'pyx_battery {s.k} {s.delta} "{s.fmt}"' for s in campaign.load()],
        f'mata: pyx_export("{(OUT / "twin_battery.csv").as_posix()}", R)',
    ]
    for fmt in DISPLAY_FORMATS:  # the same cells, written with a display format and read back
        name = f"twin_display_{fmt.strip('%').replace('.', '_')}.csv"
        lines.append(
            f'mata: pyx_export("{(OUT / name).as_posix()}", strtoreal(strofreal(R, "{fmt}")))'
        )
    lines.append("exit, clear")
    path = OUT / "self_agreement.do"
    path.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    return path


def run_stata(do_file: Path) -> None:
    """Run the do-file with Stata 17 in batch, from PowerShell only."""
    runner = HERE / "run_stata.ps1"
    subprocess.run(
        ["powershell", "-NoProfile", "-File", str(runner), "-Do", str(do_file), "-Quiet"],
        check=True,
        timeout=1800,
    )
    log = do_file.with_name(f"_run_{do_file.stem}.log").read_text(
        encoding="utf-8", errors="replace"
    )
    failing = [ln for ln in log.splitlines() if ln.startswith("r(")]
    if failing:
        raise RuntimeError(f"Stata reported errors: {failing[:3]}")


def summarize(errors: np.ndarray, columns: list[int] | None = None) -> dict[str, Any]:
    """Maxima of a ``(scenarios, quantities)`` error matrix, overall and for the worst columns."""
    keep = list(range(errors.shape[1])) if columns is None else columns
    errors = errors[:, keep]
    per_column = np.max(errors, axis=0)
    order = np.argsort(per_column)[::-1][:5]
    return {
        "max": float(per_column.max()),
        "median_of_column_maxima": float(np.median(per_column)),
        "cells_above_1e-12": int((errors > 1e-12).sum()),
        "cells": int(errors.size),
        "worst_columns": {campaign.COLUMNS[keep[i]]: float(per_column[i]) for i in order},
    }


def measure() -> dict[str, Any]:
    """Compare the three pairs of matrices (needs ``validation/out`` filled by Stata)."""
    golden = np.array(
        json.loads(GOLDEN.read_text(encoding="utf-8"))["cases"]["campaign"]["battery"], dtype=float
    )
    twin_run = read_matrix(OUT / "twin_battery.csv")
    package = np.vstack([campaign.python_row(s.frame, s.delta, s.fmt) for s in campaign.load()])
    result: dict[str, Any] = {
        "frozen_on": date.today().isoformat(),
        "scenarios": int(golden.shape[0]),
        "quantities": int(golden.shape[1]),
        "error": "|a - b| / max(|b|, 1)",
        "stata_vs_twin": summarize(error(twin_run, golden), UNIT_ORDER_FREE),
        "package_vs_stata": summarize(error(package, golden)),
        "display": {},
    }
    for fmt in DISPLAY_FORMATS:
        name = f"twin_display_{fmt.strip('%').replace('.', '_')}.csv"
        shown = read_matrix(OUT / name)
        result["display"][fmt] = summarize(error(shown, twin_run))
    return result


def main(argv: list[str]) -> int:
    """Run Stata and freeze the result, or (``--check``) read the frozen result."""
    if "--check" not in argv:
        run_stata(write_do())
        result = measure()
        TARGET.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8", newline="\n")
        print(f"wrote {TARGET}")
    else:
        result = json.loads(TARGET.read_text(encoding="utf-8"))
    for name in ("stata_vs_twin", "package_vs_stata"):
        print(f"{name:<18} max error {result[name]['max']:.2e}")
    for fmt, entry in result["display"].items():
        print(f"display {fmt:<9} max error {entry['max']:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
