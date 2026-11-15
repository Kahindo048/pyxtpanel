"""Run validation cases in Stata and freeze the reference results.

    python -m validation.run xtset            # write data + do-file, run Stata, freeze
    python -m validation.run xtset --no-stata # only (re)write the do-file

The do-file is generated from the registered cases (``validation/registry.py``).  Stata is run
in batch **from PowerShell only** by ``run_stata.ps1``.  Frozen results go to
``tests/golden/<domain>_golden.json`` with the Stata version, so the tests run without Stata.
"""

from __future__ import annotations

import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from validation import panels
from validation.compare import load_stata, read_matrix
from validation.registry import MATA_EXPORT, PRELUDES, Case, of_domain

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STATA_DIR, OUT_DIR, DATA_DIR = HERE / "stata", HERE / "out", HERE / "data"
GOLDEN_DIR = ROOT / "tests" / "golden"


def golden_path(domain: str) -> Path:
    """Frozen results of one domain."""
    return GOLDEN_DIR / f"{domain}_golden.json"


def ensure_dataset(name: str) -> None:
    """Write ``validation/data/<name>.dta`` for the simulated and hostile panels."""
    if (DATA_DIR / f"{name}.dta").exists() and not name.startswith("sim_"):
        return
    key = name.removeprefix("sim_")
    if key in panels.SPECS:
        panels.write_dta(name.removeprefix("sim_"))
    elif key in panels.HOSTILE_REFUSED:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        panels.HOSTILE_REFUSED[key].to_stata(DATA_DIR / f"{name}.dta", write_index=False)
    elif not (DATA_DIR / f"{name}.dta").exists():
        raise FileNotFoundError(f"unknown dataset {name!r}")


def render_do(domain: str) -> str:
    """Do-file that runs every case of ``domain`` and exports the objects."""
    lines = [
        "version 17",
        "set more off",
        "clear all",
        "",
        MATA_EXPORT,
        'mata: pyx_export("../out/_meta__stata.csv", (st_numscalar("c(stata_version)") \\ st_numscalar("c(MP)")))',
        "",
    ]
    if domain in PRELUDES:
        lines += [PRELUDES[domain], ""]
    for case in of_domain(domain):
        lines += [
            f"* ---- case {case.tag}",
            f'use "../data/{case.dataset}.dta", clear',
            *case.stata,
            "",
        ]
    lines.append("exit, clear")
    return "\n".join(lines) + "\n"


def write_do(domain: str) -> Path:
    """Write ``validation/stata/<domain>_cases.do`` (and the datasets it needs)."""
    STATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for case in of_domain(domain):
        ensure_dataset(case.dataset)
    path = STATA_DIR / f"{domain}_cases.do"
    path.write_text(render_do(domain), encoding="ascii", newline="\n")
    return path


def run_stata(do_file: Path) -> None:
    """Run a do-file with Stata 17 MP in batch, from PowerShell."""
    runner = HERE / "run_stata.ps1"
    subprocess.run(
        ["powershell", "-NoProfile", "-File", str(runner), "-Do", str(do_file), "-Quiet"],
        check=True,
        timeout=1800,
    )
    log = do_file.with_name(f"_run_{do_file.stem}.log").read_text(
        encoding="utf-8", errors="replace"
    )
    if "r(" in log and "\nr(" in log:
        failing = [ln for ln in log.splitlines() if ln.startswith("r(")]
        raise RuntimeError(f"Stata reported errors: {failing[:3]}")


def freeze(domain: str) -> Path:
    """Write ``tests/golden/<domain>_golden.json`` from the exports in ``validation/out``."""
    meta = read_matrix(OUT_DIR / "_meta__stata.csv").ravel()
    master = {
        "meta": {
            "stata": f"Stata {int(meta[0])} {'MP' if meta[1] else 'SE/BE'} (batch, PowerShell)",
            "frozen_on": dt.date.today().isoformat(),
            "domain": domain,
        },
        "cases": {
            case.tag: {k: v.tolist() for k, v in load_stata(OUT_DIR, case.tag).items()}
            for case in of_domain(domain)
        },
    }
    path = golden_path(domain)
    path.write_text(json.dumps(master, indent=1), encoding="utf-8", newline="\n")
    return path


def python_side(case: Case) -> dict[str, np.ndarray]:
    """Objects computed by the package for one case."""
    return {k: np.asarray(v, dtype=float) for k, v in case.python().items()}


def main(argv: list[str]) -> None:
    domain = argv[0]
    panels.freeze_panels() if "--refreeze-panels" in argv else None
    do_file = write_do(domain)
    print(f"wrote {do_file}")
    if "--no-stata" not in argv:
        run_stata(do_file)
        print(f"froze {freeze(domain)}")


if __name__ == "__main__":
    main(sys.argv[1:])
