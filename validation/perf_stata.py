"""Paired speed measure against Stata 17 MP (I.3.1).

    python -m validation.perf_stata                  # 100 000 and 1 000 000 rows
    python -m validation.perf_stata --rows 1000000   # one size

The same simulated panel (``tools.perf.simulate``) is written to a ``.dta``; a do-file times each
operation with Stata's own ``timer`` and writes the seconds to a CSV; the package then times the
equivalent operation on the same frame, best of ``--repeat``.  Stata runs through PowerShell, as
the rest of the validation.  The result is printed as a Markdown table; ``--write`` saves it.

Stata operations: ``xtset id t``; ``gen double lx = L.x``; ``gen double dx = D.x``; the within
deviation by ``bysort id: egen mean()``; ``xtsum``; ``tsfill``; ``xtdescribe``.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd
from tools.perf import operations, simulate, timed

import pyxtpanel as xt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "out"

#: Stata statements and the package operation measured against them.
PAIRS: list[tuple[str, list[str], str]] = [
    ("xtset", ["xtset id t"], "xtset"),
    ("lag", ["gen double lx = L.x"], "lag"),
    ("diff", ["gen double dx = D.x"], "diff"),
    (
        "within (unit)",
        [
            "bysort id: egen double mx = mean(x)",
            "bysort id: egen double my = mean(y)",
            "gen double wx = x - mx",
            "gen double wy = y - my",
        ],
        "within (unit)",
    ),
    (
        "between",
        ["bysort id: egen double bx = mean(x)", "bysort id: egen double by_ = mean(y)"],
        "between",
    ),
    ("xtsum", ["quietly xtsum x y"], "xtsum"),
    ("xtdescribe", ["quietly xtdescribe"], "xtdescribe"),
    ("tsfill", ["tsfill"], "tsfill"),
]


def stata_do(data: Path, result: Path) -> str:
    """The do-file: each block is timed on a fresh copy of the declared data where it needs one."""
    lines = [
        "version 17",
        "set more off",
        f'use "{data}", clear',
        'file open fh using "' + str(result) + '", write replace text',
        'file write fh "operation,seconds" _n',
    ]
    for number, (name, statements, _) in enumerate(PAIRS, start=1):
        # xtset needs the shuffled data: reload; the others run on the declared panel
        if name == "xtset":
            lines += [f'use "{data}", clear']
        elif name == "lag":
            lines += [f'use "{data}", clear', "quietly xtset id t"]
        lines += [f"timer clear {number}", f"timer on {number}"]
        lines += statements
        lines += [
            f"timer off {number}",
            f"quietly timer list {number}",
            f'file write fh "{name}," %12.6f (r(t{number})) _n',
        ]
    lines += ["file close fh", "exit, clear"]
    return "\r\n".join(lines) + "\r\n"


def run_stata(rows: int) -> dict[str, float]:
    OUT.mkdir(exist_ok=True)
    frame = simulate(rows)
    data = OUT / f"perf_{rows}.dta"
    frame.to_stata(data, write_index=False, version=118)
    result = OUT / f"perf_{rows}_stata.csv"
    do = OUT / f"perf_{rows}.do"
    do.write_text(stata_do(data, result), encoding="ascii")
    subprocess.run(
        [
            "powershell",
            "-NoProfile",
            "-File",
            str(ROOT / "run_stata.ps1"),
            "-Do",
            str(do),
            "-Quiet",
        ],
        check=True,
    )
    table = pd.read_csv(result)
    return dict(zip(table["operation"], table["seconds"], strict=True))


def run_package(rows: int, repeat: int) -> dict[str, float]:
    frame = simulate(rows)
    ops = operations()
    out = {"xtset": timed(lambda: xt.xtset(frame, "id", "t"), repeat)}
    panel = xt.xtset(frame, "id", "t")
    for _, _, key in PAIRS:
        if key != "xtset":
            ops[key](panel)
            out[key] = timed(lambda k=key: ops[k](panel), repeat)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--rows", type=int, nargs="+", default=[100_000, 1_000_000])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    lines = [
        "| Operation | Rows | Stata 17 MP (s) | pyxtpanel (s) | Ratio Stata / pyxtpanel |",
        "|---|---:|---:|---:|---:|",
    ]
    for rows in args.rows:
        stata = run_stata(rows)
        package = run_package(rows, args.repeat)
        for name, _, key in PAIRS:
            s, p = stata[name], package[key]
            lines.append(f"| {name} | {rows:,} | {s:.3f} | {p:.3f} | {s / p:.1f} |")
    text = "\n".join(lines)
    print(text)
    if args.write:
        args.write.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
