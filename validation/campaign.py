"""Random scenarios run in Stata and in the package (J.4.1).

The fixed panels of ``validation/panels.py`` are chosen by hand.  A campaign draws many more, at
random, in the dimensions where a hand-written list is thin: size (one unit, one period, a few
rows), holes, entry and exit, units with a single row, missing values (from none to all), constant
variables, a repeated (unit, period) pair, the time variable (integers with any origin and step,
Stata weeks, months and quarters), the identifier (any integers, in any order) and the order of the rows.

Every scenario is generated from the seed ``(SEED, k)``, then **frozen** with the others in
``tests/golden/data/campaign.csv.gz`` (and ``campaign_meta.json``): the tests never depend on a
random generator.  Stata runs one battery of commands on each scenario
(``validation/cases/campaign.py``); the package runs the same battery
(:func:`python_row`); the two matrices must agree.  The same frozen scenarios also feed
``tests/properties/test_campaign_invariance.py``, which relabels the identifiers, shifts the time
origin and shuffles the rows without Stata.

    python -m validation.campaign          # regenerate and freeze the scenarios (rarely)
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "tests" / "golden" / "data"
DATA_FILE = FROZEN / "campaign.csv.gz"
META_FILE = FROZEN / "campaign_meta.json"
DTA = Path(__file__).resolve().parent / "data" / "sim_campaign.dta"

SEED = 20261118
N_SCENARIOS = 160

#: Stata quarter and month numbers of a start date, in a range that is far from the epoch.
FORMATS = {"": 0.48, "%td": 0.12, "%tw": 0.07, "%tm": 0.12, "%tq": 0.11, "%th": 0.05, "%ty": 0.05}
STEPS = (1, 1, 1, 1, 2, 3, 5, 7)


@dataclass(frozen=True)
class Scenario:
    """One panel and how it is declared (``fmt`` is a Stata format for Stata-coded times)."""

    k: int
    frame: pd.DataFrame  # columns id, t, x, z
    delta: int
    fmt: str
    kind: str  # what was drawn, for the failure messages


# ---------------------------------------------------------------- drawing


def _size(rng: np.random.Generator, choices: list[tuple[int, int, float]]) -> int:
    weights = np.array([c[2] for c in choices])
    low, high, _ = choices[int(rng.choice(len(choices), p=weights / weights.sum()))]
    return int(rng.integers(low, high + 1))


def _ranks(rng: np.random.Generator, n_units: int, n_periods: int) -> tuple[list[np.ndarray], str]:
    """The calendar ranks each unit is observed in, and the name of the pattern."""
    pattern = str(rng.choice(["balanced", "holes", "entry_exit", "mixed", "ragged", "singletons"]))
    out = []
    for _ in range(n_units):
        ranks = np.arange(n_periods)
        if pattern == "holes":
            ranks = ranks[rng.random(n_periods) > float(rng.choice([0.05, 0.2, 0.5]))]
        elif pattern == "entry_exit":
            a = int(rng.integers(0, max(n_periods // 2, 1)))
            b = int(rng.integers(max(n_periods // 2, 1), n_periods + 1))
            ranks = ranks[a:b]
        elif pattern == "mixed":
            a = int(rng.integers(0, max(n_periods // 3, 1)))
            ranks = ranks[a:][rng.random(n_periods - a) > 0.15]
        elif pattern == "ragged":
            ranks = ranks[rng.random(n_periods) < float(rng.uniform(0.1, 0.9))]
        elif pattern == "singletons" and rng.random() < 0.4:
            ranks = rng.choice(n_periods, size=1)
        out.append(np.sort(ranks))
    if sum(len(r) for r in out) == 0:
        out[0] = np.array([int(rng.integers(0, n_periods))])
    return out, pattern


def _values(rng: np.random.Generator, n: int) -> tuple[np.ndarray, str]:
    """One variable: normal, constant, mostly or entirely missing, large or tiny scale."""
    style = str(
        rng.choice(
            ["normal", "normal", "normal", "missing", "constant", "all_missing", "scaled"],
            p=[0.3, 0.2, 0.2, 0.15, 0.05, 0.02, 0.08],
        )
    )
    if style == "constant":
        values = np.full(n, round(float(rng.normal()), 6))
    else:
        scale = float(rng.choice([1e-3, 1e6])) if style == "scaled" else 1.0
        values = (rng.normal(size=n) * scale).round(6)
    if style == "all_missing":
        values[:] = np.nan
    elif style == "missing":
        values[rng.random(n) < float(rng.choice([0.05, 0.3, 0.7]))] = np.nan
    elif rng.random() < 0.3:
        values[rng.random(n) < 0.08] = np.nan
    return values, style


def draw(k: int) -> Scenario:
    """Scenario ``k`` of the campaign."""
    rng = np.random.default_rng([SEED, k])
    n_units = _size(rng, [(1, 1, 0.05), (2, 5, 0.25), (6, 30, 0.55), (31, 60, 0.15)])
    n_periods = _size(rng, [(1, 1, 0.05), (2, 5, 0.25), (6, 25, 0.55), (26, 40, 0.15)])
    ranks, pattern = _ranks(rng, n_units, n_periods)
    labels = rng.choice(np.arange(-50, 5000), size=n_units, replace=False)
    unit = np.concatenate([np.full(len(r), labels[i]) for i, r in enumerate(ranks)])
    rank = np.concatenate(ranks)

    fmt = str(rng.choice(list(FORMATS), p=list(FORMATS.values())))
    delta = int(rng.choice(STEPS))
    if fmt:  # Stata date number (day, month, quarter, half-year, year): any start, any step
        time = int(rng.integers(-40, 300)) + delta * rank
    else:
        time = int(rng.integers(-30, 2000)) + delta * rank
    x, x_style = _values(rng, len(unit))
    z, z_style = _values(rng, len(unit))
    frame = pd.DataFrame({"id": unit, "t": time, "x": x, "z": z})
    kind = f"{n_units}x{n_periods} {pattern}, x {x_style}, z {z_style}, t {fmt} delta {delta}"
    if rng.random() < 0.03:  # a repeated (unit, period) pair: Stata refuses with 451
        frame = pd.concat([frame, frame.sample(1, random_state=int(rng.integers(1 << 30)))])
        kind += ", repeated pair"
    if rng.random() < 0.5:
        frame = frame.iloc[rng.permutation(len(frame))]
        kind += ", shuffled"
    return Scenario(k, frame.reset_index(drop=True), delta, fmt, kind)


# ---------------------------------------------------------------- freezing


def freeze() -> None:
    """Draw every scenario and write the frozen files."""
    FROZEN.mkdir(parents=True, exist_ok=True)
    scenarios = [draw(k) for k in range(N_SCENARIOS)]
    stacked = pd.concat([s.frame.assign(scenario=s.k) for s in scenarios], ignore_index=True)[
        ["scenario", "id", "t", "x", "z"]
    ]
    stacked.to_csv(DATA_FILE, index=False, compression="gzip")
    meta = {str(s.k): {"delta": s.delta, "fmt": s.fmt, "kind": s.kind} for s in scenarios}
    META_FILE.write_text(json.dumps(meta, indent=1), encoding="utf-8", newline="\n")


def load() -> list[Scenario]:
    """The frozen scenarios."""
    stacked = pd.read_csv(DATA_FILE)
    meta = json.loads(META_FILE.read_text(encoding="utf-8"))
    return [
        Scenario(
            int(k),
            stacked[stacked["scenario"] == int(k)][["id", "t", "x", "z"]].reset_index(drop=True),
            int(info["delta"]),
            str(info["fmt"]),
            str(info["kind"]),
        )
        for k, info in meta.items()
    ]


def write_dta() -> Path:
    """``validation/data/sim_campaign.dta``: every scenario, one after the other."""
    DTA.parent.mkdir(parents=True, exist_ok=True)
    stacked = pd.read_csv(DATA_FILE)
    stacked["scenario"] = stacked["scenario"].astype("int32")
    stacked["id"] = stacked["id"].astype("int32")
    stacked["t"] = stacked["t"].astype("int32")
    stacked.to_stata(DTA, write_index=False, version=118)
    return DTA


# ---------------------------------------------------------------- the battery

CORR_KINDS = ("overall", "within", "between")
XTSUM = [
    "N",
    "n",
    "Tbar",
    "mean",
    "sd",
    "min",
    "max",
    "sd_b",
    "min_b",
    "max_b",
    "sd_w",
    "min_w",
    "max_w",
]
OPERATORS = {"L": "L.x", "L2": "L2.x", "F": "F.x", "D": "D.x", "D2": "D2.x", "S2": "S2.x"}

#: Names of the columns of the battery, in order; the Stata program fills the same ones.
COLUMNS: list[str] = [
    "rc", "tmin", "tmax", "tdelta", "gaps", "bal", "nobs", "nunits", "nperiods",
    *[f"xsum_{n}" for n in XTSUM],
    *[f"zsum_{n}" for n in XTSUM],
    *[f"{op}_{stat}" for op in OPERATORS for stat in ("n", "sum", "wsum")],
    "within_ssq", "within_wsum", "between_sum", "between_wsum",
    "desc_rc", "desc_units", "desc_tmin", "desc_tmax",
    "fill_nobs", "fill_nx", "fillfull_nobs",
    "trans_n", "trans_rows", "trans_cols",
    *[f"corr_{kind}_{what}" for kind in CORR_KINDS for what in ("r", "n")],
]  # fmt: skip

#: Columns that do not change when identifiers are relabelled, time shifted or rows shuffled
#: (the others hold a value of the time variable).
LABEL_FREE = [c for c in COLUMNS if c not in ("tmin", "tmax")]


def _triple(values: np.ndarray) -> list[float]:
    """Count, sum and sum of position x value of the observed values (positions from 1)."""
    keep = ~np.isnan(values)
    position = np.arange(1, len(values) + 1, dtype=float)
    return [
        float(keep.sum()),
        float(values[keep].sum()),
        float((values[keep] * position[keep]).sum()),
    ]


def _fill(row: np.ndarray, values: dict[str, float]) -> None:
    for name, value in values.items():
        row[COLUMNS.index(name)] = value


def _structure(row: np.ndarray, panel) -> None:  # type: ignore[no-untyped-def]
    first, last = int(panel.index.rank.min()), int(panel.index.rank.max())
    _fill(
        row,
        {
            "tmin": float(panel.axis.stata_values(first)[0]),
            "tmax": float(panel.axis.stata_values(last)[0]),
            "tdelta": float(panel.delta),
            "gaps": float(panel.has_gaps),
            "bal": {"strongly balanced": 1, "weakly balanced": 2, "unbalanced": 3}[panel.balanced],
            "nobs": panel.n_obs,
            "nunits": panel.n_units,
            "nperiods": panel.n_periods,
        },
    )


def _summaries(row: np.ndarray, panel) -> None:  # type: ignore[no-untyped-def]
    from pyxtpanel import PyxtpanelError

    for name in ("x", "z"):
        try:
            scalars = panel.xtsum(name).scalars.iloc[0]
        except PyxtpanelError:  # noqa: PERF203 - Stata: r(2000), no observation
            continue
        _fill(row, {f"{name}sum_{stat}": float(scalars[stat]) for stat in XTSUM})


def _transformations(row: np.ndarray, panel) -> None:  # type: ignore[no-untyped-def]
    order = panel.index.source_row
    ops = panel.operate(*OPERATORS.values()).iloc[order]
    for op, term in OPERATORS.items():
        triple = _triple(ops[term].to_numpy())
        _fill(
            row, {f"{op}_{stat}": v for stat, v in zip(("n", "sum", "wsum"), triple, strict=True)}
        )
    within = panel.within("x").iloc[order].to_numpy()
    between = panel.between("x").iloc[order].to_numpy()
    position = np.arange(1, len(within) + 1, dtype=float)
    seen, kept = ~np.isnan(within), ~np.isnan(between)
    _fill(
        row,
        {
            "within_ssq": float((within[seen] ** 2).sum()),
            "within_wsum": float((within[seen] * position[seen]).sum()),
            "between_sum": float(between[kept].sum()),
            "between_wsum": float((between[kept] * position[kept]).sum()),
        },
    )


def _description(row: np.ndarray, panel) -> None:  # type: ignore[no-untyped-def]
    from pyxtpanel import DataError

    try:
        panel.xtdescribe()
    except DataError as error:  # Stata: r(459), a panel observed in a single period
        _fill(row, {"desc_rc": float(error.stata_code or np.nan)})
        return
    counts = panel.counts
    _fill(
        row,
        {
            "desc_rc": 0.0,
            "desc_units": float(panel.n_units),
            "desc_tmin": float(np.min(counts)),
            "desc_tmax": float(np.max(counts)),
        },
    )


def _tables_and_correlations(row: np.ndarray, panel, frame, delta, fmt) -> None:  # type: ignore[no-untyped-def]
    from pyxtpanel import PyxtpanelError, xtset

    z = frame["z"].to_numpy()
    discrete = np.select([z < -0.5, z < 0.5, z >= 0.5], [0.0, 1.0, 2.0], default=np.nan)
    try:
        coded = xtset(frame.assign(dz=discrete), "id", "t", delta=delta, format=fmt or None)
        trans = coded.xttrans("dz")
        _fill(
            row,
            {"trans_n": trans.n_trans, "trans_rows": trans.n_rows, "trans_cols": trans.n_cols},
        )
    except PyxtpanelError:
        pass  # Stata: nothing to tabulate
    try:
        corr = panel.xtcorr(["x", "z"])
    except PyxtpanelError:
        return
    for kind in CORR_KINDS:
        n_points, r = float(corr.n[kind].loc["x", "z"]), float(corr.corr[kind].loc["x", "z"])
        # DIVERGENCE: Stata's correlate refuses fewer than 2 observations (r(2001)) and stores
        # nothing; the package reports the number of points with a missing correlation
        # (tests/unit/test_campaign_findings.py).  A variable constant up to rounding has a
        # missing correlation here, and the Stata side does not read the noise it leaves.
        if n_points >= 2 and not np.isnan(r):
            _fill(row, {f"corr_{kind}_r": r, f"corr_{kind}_n": n_points})


def python_row(frame: pd.DataFrame, delta: int, fmt: str) -> np.ndarray:
    """The battery as the package computes it (``NaN`` where Stata would refuse or have no value)."""
    from pyxtpanel import DeclarationError, xtset

    row = np.full(len(COLUMNS), np.nan)
    try:
        panel = xtset(frame, "id", "t", delta=delta, format=fmt or None)
    except DeclarationError as error:
        row[0] = error.stata_code or np.nan
        return row
    row[0] = 0
    _structure(row, panel)
    _summaries(row, panel)
    _transformations(row, panel)
    _description(row, panel)
    _tables_and_correlations(row, panel, frame, delta, fmt)
    filled = panel.tsfill()
    _fill(
        row,
        {
            "fill_nobs": float(filled.n_obs),
            "fill_nx": float(filled.to_frame()["x"].notna().sum()),
            "fillfull_nobs": float(panel.tsfill(full=True).n_obs),
        },
    )
    return row


# ---------------------------------------------------------------- the Stata side


def _corr_block(kind: str, put) -> list[str]:  # type: ignore[no-untyped-def]
    """Stata lines: the correlation of x and z on the rows where both are observed."""
    lines = [
        "        quietly use `base', clear",
        "        quietly xtset id t, delta(`delta')",
        "        quietly keep if !missing(x, z)",
    ]
    if kind != "overall":
        option = "fe" if kind == "within" else "be"
        lines += [
            f"        capture quietly xtdata x z, {option} clear",
            "        local rcdata = _rc",
        ]
    else:
        lines += ["        local rcdata = 0"]
    lines += [
        "        // a variable constant up to rounding has no correlation: Stata's xtdata leaves noise",
        "        // of 1e-16 in it and correlate reports a number made of that noise, so the battery",
        "        // does not read it (the package returns a missing value, see describe/corr.py)",
        "        if `rcdata' == 0 {",
        "            quietly summarize x",
        "            local degenerate = missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))",
        "            quietly summarize z",
        "            local degenerate = `degenerate' | missing(r(sd)) | r(sd) <= 1e-12 * max(abs(r(min)), abs(r(max)))",
        "            capture quietly correlate x z",
        "            if _rc == 0 & !`degenerate' {",
        "    " + put(f"corr_{kind}_r", "r(rho)"),
        "    " + put(f"corr_{kind}_n", "r(N)"),
        "            }",
        "        }",
    ]
    return lines


def stata_program() -> str:
    """The Stata program that fills one row of the battery for the scenario it is given."""
    at = {name: i + 1 for i, name in enumerate(COLUMNS)}

    def put(name: str, expression: str) -> str:
        return f"    matrix row[1, {at[name]}] = {expression}"

    lines = [
        "capture program drop pyx_battery",
        "program define pyx_battery",
        "    args k delta fmt",
        "    preserve",
        "    quietly keep if scenario == `k'",
        f"    matrix row = J(1, {len(COLUMNS)}, .)",
        '    if "`fmt\'" != "" format t `fmt\'',
        "    capture xtset id t, delta(`delta')",
        put("rc", "_rc"),
        "    if _rc == 0 {",
        put("tmin", "r(tmin)"),
        put("tmax", "r(tmax)"),
        put("tdelta", "r(tdelta)"),
        put("gaps", "r(gaps)"),
        put(
            "bal",
            'cond("`r(balanced)\'" == "strongly balanced", 1, '
            'cond("`r(balanced)\'" == "weakly balanced", 2, 3))',
        ),
        "        scalar nobs_ = _N",
        "        tempvar f1 f2",
        "        bysort id: gen byte `f1' = (_n == 1)",
        "        quietly summarize `f1'",
        "        scalar nunits_ = r(sum)",
        "        bysort t: gen byte `f2' = (_n == 1)",
        "        quietly summarize `f2'",
        "        scalar nper_ = r(sum)",
        "        sort id t",
        put("nobs", "nobs_"),
        put("nunits", "nunits_"),
        put("nperiods", "nper_"),
    ]
    for name in ("x", "z"):
        lines.append(f"        capture quietly xtsum {name}")
        lines.append("        if _rc == 0 {")
        lines.extend("    " + put(f"{name}sum_{stat}", f"r({stat})") for stat in XTSUM)
        lines.append("        }")
    for op, term in OPERATORS.items():
        lines += [
            f"        tempvar o_{op}",
            f"        quietly gen double `o_{op}' = {term}",
            f'        mata: st_matrix("t3", pyx_triple("`o_{op}\'"))',
            *[
                "    " + put(f"{op}_{stat}", f"t3[1, {j + 1}]")
                for j, stat in enumerate(("n", "sum", "wsum"))
            ],
        ]
    lines += [
        "        tempvar mx dev bx",
        "        quietly bysort id: egen double `mx' = mean(x)",
        "        quietly gen double `dev' = x - `mx'",
        "        quietly gen double `bx' = `mx' if !missing(x)",
        '        mata: st_matrix("w2", pyx_pair("`dev\'"))',
        '        mata: st_matrix("b2", pyx_pair("`bx\'"))',
        "    " + put("within_ssq", "w2[1, 1]"),
        "    " + put("within_wsum", "w2[1, 2]"),
        "    " + put("between_sum", "b2[1, 3]"),
        "    " + put("between_wsum", "b2[1, 2]"),
        "        capture quietly xtdescribe",
        "    " + put("desc_rc", "_rc"),
        "        if _rc == 0 {",
        "    " + put("desc_units", "r(N)"),
        "    " + put("desc_tmin", "r(min)"),
        "    " + put("desc_tmax", "r(max)"),
        "        }",
        "        tempfile base",
        "        quietly save `base'",
        "        capture tsfill",
        "        if _rc == 0 {",
        "    " + put("fill_nobs", "_N"),
        "            quietly count if !missing(x)",
        "    " + put("fill_nx", "r(N)"),
        "        }",
        "        quietly use `base', clear",
        "        quietly xtset id t, delta(`delta')",
        "        capture tsfill, full",
        "        if _rc == 0 {",
        "    " + put("fillfull_nobs", "_N"),
        "        }",
        "        quietly use `base', clear",
        "        quietly xtset id t, delta(`delta')",
        "        quietly generate byte dz = 0 if z < -0.5",
        "        quietly replace dz = 1 if z >= -0.5 & z < 0.5",
        "        quietly replace dz = 2 if z >= 0.5 & !missing(z)",
        "        capture quietly xttrans dz",
        "        if _rc == 0 {",
        "    " + put("trans_n", "r(n_trans)"),
        "    " + put("trans_rows", "r(n_rows)"),
        "    " + put("trans_cols", "r(n_cols)"),
        "        }",
        *[line for kind in CORR_KINDS for line in _corr_block(kind, put)],
        "    }",
        '    mata: R = R \\ st_matrix("row")',
        "    restore",
        "end",
    ]
    return "\n".join(lines)


#: Mata helpers: (count, sum, sum of position x value) of the observed values, and the pair
#: (sum of squares, sum of position x value, sum) used for the within and between statistics.
MATA = """\
mata:
real rowvector pyx_triple(string scalar v)
{
    real colvector x, k, keep
    x = st_data(., v); k = (1::rows(x)); keep = (x :< .)
    return((sum(keep), sum(select(x, keep)), sum(select(x :* k, keep))))
}
real rowvector pyx_pair(string scalar v)
{
    real colvector x, k, keep
    x = st_data(., v); k = (1::rows(x)); keep = (x :< .)
    return((sum(select(x :^ 2, keep)), sum(select(x :* k, keep)), sum(select(x, keep))))
}
end
"""


if __name__ == "__main__":
    freeze()
    print(f"froze {N_SCENARIOS} scenarios in {DATA_FILE}")
