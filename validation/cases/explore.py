"""Variation, spells and event time against hand-written Stata logic (K.7.1).

Stata has no command for these tools, so each is checked against the counts that ``egen``,
``bysort`` and lags give on the same simulated panels.  The discrete variable ``d`` is
``x > 0``.
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "explore"
PANELS = (
    "balanced",
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "unsorted_ids",
    "one_obs_unit",
    "re_effects",
)
VARIABLES = ("x", "z")


def _option(name: str) -> str:
    delta = panels.SPECS[name].delta
    return f", delta({delta})" if delta != 1 else ""


def _panel(name: str):  # type: ignore[no-untyped-def]
    from pyxtpanel import xtset

    frame = panels.load(name).copy()
    frame["d"] = np.where(frame["x"].isna(), np.nan, (frame["x"] > 0).astype(float))
    return xtset(frame, "id", "t", delta=panels.SPECS[name].delta)


# ---------------------------------------------------------------- variation


def _variation(name: str) -> None:
    tag = f"variation_{name}"
    lines: list[str] = []
    scalars: list[str] = []
    for v in VARIABLES:
        lines += [
            "preserve",
            f"quietly keep if !missing({v})",
            f"bysort id: egen double mn = min({v})",
            f"bysort id: egen double mx = max({v})",
            "bysort id: gen byte u1 = (_n == 1)",
            "bysort id: gen long cnt = _N",
            "quietly count if u1 & mx > mn",
            f"scalar uv_{v} = r(N)",
            "quietly count if u1 & cnt == 1",
            f"scalar us_{v} = r(N)",
            f"bysort t: egen double tmn = min({v})",
            f"bysort t: egen double tmx = max({v})",
            "bysort t: gen byte p1 = (_n == 1)",
            "quietly count if p1 & tmx > tmn",
            f"scalar pv_{v} = r(N)",
            f"bysort {v}: gen byte d1 = (_n == 1)",
            "quietly count if d1",
            f"scalar nd_{v} = r(N)",
            "restore",
        ]
        scalars += [f"uv_{v}", f"us_{v}", f"pv_{v}", f"nd_{v}"]
    lines.append(export_scalars(tag, "counts", scalars))

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        table = _panel(name).variation(list(VARIABLES)).table
        out = []
        for v in VARIABLES:
            row = table.loc[v]
            out += [
                row["units_varying"],
                row["units_single"],
                row["periods_varying"],
                row["n_distinct"],
            ]
        return {"counts": np.array(out, dtype=float).reshape(-1, 1)}


# ---------------------------------------------------------------- spells


def _spells(name: str) -> None:
    tag = f"spells_{name}"
    lines = [
        f"xtset id t{_option(name)}",
        "generate byte d = cond(missing(x), ., x > 0)",
        "generate byte ld = L.d",
        "generate byte start = !missing(d) & (missing(ld) | ld != d)",
        "generate byte change = !missing(d) & !missing(ld) & ld != d",
        "quietly count if start",
        "scalar n_spells = r(N)",
        "quietly count if change",
        "scalar n_changes = r(N)",
        "bysort id: egen byte any_change = max(change)",
        "bysort id: gen byte u1 = (_n == 1)",
        "quietly count if u1 & any_change == 1",
        "scalar n_units_changing = r(N)",
        export_scalars(tag, "counts", ["n_spells", "n_changes", "n_units_changing"]),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        result = _panel(name).spells("d")
        counts = [len(result.spells), result.n_changes, result.n_units_changing]
        return {"counts": np.array(counts, dtype=float).reshape(-1, 1)}


# ---------------------------------------------------------------- event time


def _event(name: str) -> None:
    tag = f"event_{name}"
    delta = panels.SPECS[name].delta
    lines = [
        f"xtset id t{_option(name)}",
        "generate byte d = cond(missing(x), ., x > 0)",
        "generate double on = cond(d == 1, t, .)",
        "bysort id: egen double ev = min(on)",
        f"generate double rel = (t - ev) / {delta}",
        "quietly summarize rel",
        "scalar n_rel = r(N)",
        "scalar sum_rel = r(sum)",
        "generate double rel2 = rel^2",
        "quietly summarize rel2",
        "scalar sumsq_rel = r(sum)",
        "bysort id: gen byte u1 = (_n == 1)",
        "quietly count if u1 & !missing(ev)",
        "scalar n_treated = r(N)",
        export_scalars(tag, "counts", ["n_rel", "sum_rel", "sumsq_rel", "n_treated"]),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        result = _panel(name).event_time(indicator="d")
        rel = result.relative.dropna().to_numpy(dtype=float)
        counts = [rel.size, rel.sum(), (rel**2).sum(), result.n_treated]
        return {"counts": np.array(counts, dtype=float).reshape(-1, 1)}


# ---------------------------------------------------------------- within outliers


def _outliers(name: str) -> None:
    tag = f"outliers_{name}"
    lines: list[str] = []
    scalars: list[str] = []
    for v in VARIABLES:
        lines += [
            f"bysort id: egen long n_{v} = count({v})",
            "* egen median() averages the two middle values in float precision; pctile() does not",
            f"bysort id: egen double med_{v} = pctile({v}), p(50)",
            f"generate double dev_{v} = abs({v} - med_{v})",
            f"bysort id: egen double mad_{v} = pctile(dev_{v}), p(50)",
            f"generate double z_{v} = 0.6745 * ({v} - med_{v}) / mad_{v} "
            f"if n_{v} >= 3 & mad_{v} > 0 & !missing({v})",
            f"generate double az_{v} = abs(z_{v})",
            f"quietly summarize az_{v}",
            f"scalar ns_{v} = r(N)",
            f"scalar sa_{v} = r(sum)",
            f"quietly count if az_{v} > 2 & !missing(az_{v})",
            f"scalar nf_{v} = r(N)",
        ]
        scalars += [f"ns_{v}", f"nf_{v}", f"sa_{v}"]
    lines.append(export_scalars(tag, "counts", scalars))

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines, atol=1e-9, rtol=1e-9)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        panel = _panel(name)
        out: list[float] = []
        for v in VARIABLES:
            result = panel.within_outliers(v, threshold=2.0, min_obs=3)
            out += [result.n_scored, result.n_flagged, float(np.nansum(np.abs(result.scores)))]
        return {"counts": np.array(out, dtype=float).reshape(-1, 1)}


# ---------------------------------------------------------------- wide form


def _wide(name: str) -> None:
    tag = f"wide_{name}"
    lines = [
        "keep id t x z",
        "reshape wide x z, i(id) j(t)",
        "scalar ncols = c(k)",
        "scalar nunits = _N",
        "egen double tot = rowtotal(x* z*)",
        "quietly summarize tot",
        "scalar total = r(sum)",
        "egen long nm = rownonmiss(x* z*)",
        "quietly summarize nm",
        "scalar nonmiss = r(sum)",
        'mata: M = editmissing(st_data(., "x* z*"), 0)',
        'mata: st_numscalar("colw", sum((1..cols(M)) :* colsum(M)))',
        'mata: st_numscalar("roww", sum((1::rows(M)) :* rowsum(M)))',
        export_scalars(tag, "counts", ["ncols", "nunits", "total", "nonmiss", "colw", "roww"]),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines, atol=1e-9, rtol=1e-9)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        wide = _panel(name).to_wide(["x", "z"])
        grid = wide.fillna(0.0).to_numpy()
        rows = np.arange(1, grid.shape[0] + 1)
        cols = np.arange(1, grid.shape[1] + 1)
        out = [
            wide.shape[1] + 1,  # Stata counts the identifier
            wide.shape[0],
            grid.sum(),
            int(wide.notna().to_numpy().sum()),
            float((cols * grid.sum(axis=0)).sum()),
            float((rows * grid.sum(axis=1)).sum()),
        ]
        return {"counts": np.array(out, dtype=float).reshape(-1, 1)}


# ---------------------------------------------------------------- gaps (tsreport)


def _gaps(name: str) -> None:
    tag = f"gaps_{name}"
    lines = [
        f"xtset id t{_option(name)}",
        "quietly tsreport, panel",
        "scalar n_gaps = r(N_gaps)",
        "scalar n_obs = r(N_obs)",
        export_scalars(tag, "counts", ["n_gaps", "n_obs"]),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        panel = _panel(name)
        out = [len(panel.gap_types().gaps), panel.n_obs]
        return {"counts": np.array(out, dtype=float).reshape(-1, 1)}


for _name in PANELS:
    _variation(_name)
    _spells(_name)
    _event(_name)
    _outliers(_name)
    _wide(_name)
    _gaps(_name)
