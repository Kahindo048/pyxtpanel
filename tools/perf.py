"""Speed and memory of the main operations on large panels (I.1).

    python -m tools.perf                       # 100 000 and 1 000 000 rows, one line per operation
    python -m tools.perf --rows 100000 --profile xtset within   # cProfile of chosen operations
    python -m tools.perf --pandas              # reference timings of the pandas equivalents
    python -m tools.perf --baseline tests/benchmarks/baseline.json   # reference of the speed guard

The panel is simulated (``T = 50`` periods, about 5% of the rows removed at random, three numeric
columns); the same seed gives the same panel.  Each operation is timed from a freshly declared
panel, best of ``--repeat`` runs, so that cached results cannot flatter it.
"""

from __future__ import annotations

import argparse
import cProfile
import gc
import io
import json
import platform
import pstats
import time
import tracemalloc
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd

import pyxtpanel as xt

T = 50


def simulate(rows: int, *, seed: int = 20261117) -> pd.DataFrame:
    """A panel of about ``rows`` rows: ``rows / T`` units, ``T`` periods, 5% of the rows removed."""
    rng = np.random.default_rng(seed)
    n_units = max(rows // T, 2)
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(n_units, dtype=np.int64), T),
            "t": np.tile(np.arange(T, dtype=np.int64), n_units),
        }
    )
    frame = frame[rng.random(len(frame)) > 0.05].reset_index(drop=True)
    effect = rng.normal(size=n_units)[frame["id"].to_numpy()]
    frame["x"] = effect + rng.normal(size=len(frame))
    frame["y"] = 0.5 * frame["x"] + effect + rng.normal(size=len(frame))
    frame["z"] = rng.normal(size=len(frame))
    # shuffle the rows: the worst case for a declaration that has to sort
    return frame.iloc[rng.permutation(len(frame))].reset_index(drop=True)


def operations() -> dict[str, Callable[[xt.PanelData], object]]:
    """Name -> function of a declared panel, in the order of the report."""
    return {
        "lag": lambda p: p.lag("x", 1),
        "diff": lambda p: p.diff("x"),
        "operate L2D.x": lambda p: p.operate("L2D.x"),
        "within (unit)": lambda p: p.within(["x", "y"]),
        "within (two-way)": lambda p: p.within(["x", "y"], effects="twoway"),
        "between": lambda p: p.between(["x", "y"]),
        "fod": lambda p: p.fod("x"),
        "group_mean": lambda p: p.group_mean("x"),
        "xtsum": lambda p: p.xtsum(["x", "y"]),
        "xtdescribe": lambda p: p.xtdescribe(),
        "xttrans (binned)": lambda p: p.xttrans("x", bins=5),
        "xtcorr": lambda p: p.xtcorr(["x", "y", "z"]),
        "gaps": lambda p: p.gaps(),
        "gap_types": lambda p: p.gap_types(),
        "tsfill": lambda p: p.tsfill(),
        "to_wide": lambda p: p.to_wide(["x"]),
        "arrays (contract)": lambda p: p.arrays(["y", "x"]),
    }


def pandas_equivalents() -> dict[str, Callable[[pd.DataFrame], object]]:
    """The pandas way of doing what some operations do, for the comparison of I.3.2."""

    def lag(f: pd.DataFrame) -> object:
        f = f.sort_values(["id", "t"])
        shifted = f.groupby("id")["x"].shift(1)
        ok = f.groupby("id")["t"].diff() == 1  # the calendar lag, not the previous row
        return shifted.where(ok)

    def within(f: pd.DataFrame) -> object:
        return f[["x", "y"]] - f.groupby("id")[["x", "y"]].transform("mean")

    return {
        "lag": lag,
        "within (unit)": within,
        "between": lambda f: f.groupby("id")[["x", "y"]].transform("mean"),
        "xtsum": lambda f: f[["x", "y"]].describe(),
        "tsfill": lambda f: f.set_index(["id", "t"]).reindex(
            pd.MultiIndex.from_product([f["id"].unique(), np.arange(T)], names=["id", "t"])
        ),
    }


def timed(function: Callable[[], object], repeat: int) -> float:
    best = float("inf")
    for _ in range(repeat):
        gc.collect()
        start = time.perf_counter()
        function()
        best = min(best, time.perf_counter() - start)
    return best


def peak_megabytes(function: Callable[[], object]) -> float:
    gc.collect()
    tracemalloc.start()
    function()
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    return peak / 2**20


def calibration(repeat: int = 7) -> float:
    """Seconds of a fixed numpy workload (sort, grouped sum, gather) on this machine, right now.

    The speed guard divides each operation by it, so that a slower or busier machine moves both
    and the ratio stays comparable with the one recorded in ``baseline.json``.
    """
    rng = np.random.default_rng(0)
    values = rng.random(1_000_000)
    groups = rng.integers(0, 20_000, 1_000_000)

    def kernel() -> object:
        order = np.argsort(values)
        return np.bincount(groups, weights=values, minlength=20_000), values[order].cumsum()

    return timed(kernel, repeat)


def guard_rows() -> int:
    """Size of the panel the speed guard uses."""
    return 200_000


def baseline() -> dict[str, object]:
    """Ratio operation / calibration of every operation on the guard panel (best of 7 runs)."""
    reference = calibration()
    runs = measure(guard_rows(), 7, memory=False)
    return {
        "rows": guard_rows(),
        "calibration_s": round(reference, 6),
        "ratios": {str(r["operation"]): round(float(r["s"]) / reference, 4) for r in runs},  # type: ignore[arg-type]
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
    }


def measure(rows: int, repeat: int, *, memory: bool) -> list[dict[str, object]]:
    frame = simulate(rows)
    out: list[dict[str, object]] = []
    declare = timed(lambda: xt.xtset(frame, "id", "t"), repeat)
    row: dict[str, object] = {
        "operation": "xtset (shuffled rows)",
        "rows": len(frame),
        "s": declare,
    }
    if memory:
        row["peak_mb"] = peak_megabytes(lambda: xt.xtset(frame, "id", "t"))
    out.append(row)
    for name, function in operations().items():
        panel = xt.xtset(frame, "id", "t")
        function(panel)  # warm the lazy parts that belong to the panel, not to the operation
        row = {
            "operation": name,
            "rows": len(frame),
            "s": timed(lambda f=function, p=panel: f(p), repeat),
        }
        if memory:
            row["peak_mb"] = peak_megabytes(lambda f=function, p=panel: f(p))
        out.append(row)
    return out


def measure_pandas(rows: int, repeat: int) -> dict[str, float]:
    frame = simulate(rows)
    return {
        name: timed(lambda f=function: f(frame), repeat)
        for name, function in pandas_equivalents().items()
    }


def profile(rows: int, names: list[str], top: int) -> str:
    frame = simulate(rows)
    ops = {"xtset": lambda p: xt.xtset(frame, "id", "t"), **operations()}
    text = io.StringIO()
    for name in names:
        function = ops[name]
        panel = xt.xtset(frame, "id", "t")
        profiler = cProfile.Profile()
        profiler.enable()
        function(panel)
        profiler.disable()
        text.write(f"\n=== {name} ({rows} rows) ===\n")
        stats = pstats.Stats(profiler, stream=text).sort_stats("cumulative")
        stats.print_stats(top)
    return text.getvalue()


def table(
    results: dict[int, list[dict[str, object]]], reference: dict[int, dict[str, float]]
) -> str:
    """Markdown table: seconds, peak memory and the pandas time where pandas has an equivalent."""
    sizes = sorted(results)
    head = ["Operation"]
    for n in sizes:
        head += [f"{n:,} rows (s)", "peak (MB)"]
        if reference:
            head.append("pandas (s)")
    lines = ["| " + " | ".join(head) + " |", "|---|" + "---:|" * (len(head) - 1)]
    for name in [str(r["operation"]) for r in results[sizes[0]]]:
        cells = [name]
        for n in sizes:
            row = next(r for r in results[n] if r["operation"] == name)
            cells.append(f"{float(row['s']):.3f}")  # type: ignore[arg-type]
            cells.append(f"{float(row['peak_mb']):.0f}" if "peak_mb" in row else "")  # type: ignore[arg-type]
            if reference:
                other = reference[n].get(name.replace(" (shuffled rows)", ""))
                cells.append("" if other is None else f"{other:.3f}")
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--rows", type=int, nargs="+", default=[100_000, 1_000_000])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--memory", action="store_true", help="peak allocation of each operation")
    parser.add_argument("--pandas", action="store_true", help="also time the pandas equivalents")
    parser.add_argument("--profile", nargs="+", metavar="OPERATION", help="cProfile of these")
    parser.add_argument("--top", type=int, default=12)
    parser.add_argument("--write", type=Path, help="write the Markdown report here")
    parser.add_argument("--baseline", type=Path, help="write the reference ratios of the guard")
    args = parser.parse_args()
    if args.baseline:
        args.baseline.write_text(json.dumps(baseline(), indent=2) + "\n", encoding="utf-8")
        print(f"wrote {args.baseline}")
        return
    if args.profile:
        for rows in args.rows:
            print(profile(rows, args.profile, args.top))
        return
    results = {rows: measure(rows, args.repeat, memory=args.memory) for rows in args.rows}
    reference = (
        {rows: measure_pandas(rows, args.repeat) for rows in args.rows} if args.pandas else {}
    )
    report = table(results, reference)
    for rows, runs in results.items():
        for run in runs:
            extra = f"  peak {run['peak_mb']:.0f} MB" if "peak_mb" in run else ""
            print(f"{rows:>9,}  {run['operation']:<24} {float(run['s']):8.3f} s{extra}")  # type: ignore[arg-type]
        for name, seconds in reference.get(rows, {}).items():
            print(f"{rows:>9,}  pandas {name:<17} {seconds:8.3f} s")
    if args.write:
        args.write.write_text(report + "\n", encoding="utf-8")
        print(f"wrote {args.write}")


if __name__ == "__main__":
    main()
