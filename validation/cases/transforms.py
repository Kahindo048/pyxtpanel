"""Within, between and random-effects quasi-demeaning against Stata (D.4.2, D.4.5).

``xtdata, fe`` and ``xtdata, be`` are the reference for ``within(add_mean=True, listwise=True)``
and ``between(listwise=True)``; the quasi-demeaning is checked against a transformation written
by hand in Stata from the ``σ_u`` and ``σ_e`` of ``xtreg, re`` (their value is frozen here and
the do-file checks that Stata still finds it).
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export

DOMAIN = "transforms"
XTDATA_PANELS = (
    "balanced",
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "quarterly",
    "unsorted_ids",
    "one_obs_unit",
    "single_unit",
    "single_period",
    "re_effects",
)


def _declared(name: str):  # type: ignore[no-untyped-def]
    from pyxtpanel import xtset

    spec = panels.SPECS[name]
    frame = panels.load(name)
    return frame, xtset(frame, "id", "t", delta=spec.delta)


def _option(name: str) -> str:
    delta = panels.SPECS[name].delta
    return f", delta({delta})" if delta != 1 else ""


def _within(name: str) -> dict[str, np.ndarray]:
    frame, panel = _declared(name)
    out = panel.within(["x", "z"], add_mean=True, listwise=True)
    rows = panel.index.source_row
    table = np.column_stack(
        [frame["id"].to_numpy(float)[rows], out["x"].to_numpy()[rows], out["z"].to_numpy()[rows]]
    )
    return {"data": table[~np.isnan(table).any(axis=1)]}


def _between(name: str) -> dict[str, np.ndarray]:
    frame, panel = _declared(name)
    out = panel.between(["x", "z"], listwise=True)
    rows = panel.index.source_row
    table = np.column_stack(
        [frame["id"].to_numpy(float)[rows], out["x"].to_numpy()[rows], out["z"].to_numpy()[rows]]
    )
    table = table[~np.isnan(table).any(axis=1)]
    first = np.unique(table[:, 0], return_index=True)[1]  # xtdata, be keeps one row per unit
    return {"data": table[np.sort(first)]}


def _register_xtdata(name: str) -> None:
    for mode, python in (("fe", _within), ("be", _between)):
        tag = f"transforms_{mode}_{name}"
        lines = [
            f"xtset id t{_option(name)}",
            f"xtdata x z, {mode} clear",
            "sort id",
            export(tag, "data", 'st_data(., ("id", "x", "z"))'),
        ]

        @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
        def _run(name: str = name, python=python) -> dict[str, np.ndarray]:  # type: ignore[no-untyped-def]
            return python(name)


for _name in XTDATA_PANELS:
    _register_xtdata(_name)

# ---------------------------------------------------------------- quasi-demeaning (D.4.5)

#: ``σ_u`` and ``σ_e`` of ``xtreg x z, re`` on ``re_effects`` (frozen; the do-file checks them).
SIGMA_U, SIGMA_E = 1.119347976677613, 0.7748551012200685


def _quasi() -> dict[str, np.ndarray]:
    frame, panel = _declared("re_effects")
    theta = panel.theta_re(SIGMA_U, SIGMA_E)
    out = panel.quasi_demean(["x", "z"], theta.to_numpy())
    rows = panel.index.source_row
    table = np.column_stack(
        [
            frame["id"].to_numpy(float)[rows],
            out["x"].to_numpy()[rows],
            out["z"].to_numpy()[rows],
            frame["id"].map(theta).to_numpy()[rows],
        ]
    )
    return {"data": table, "sigma_drift": np.zeros((4, 1))}


_QUASI_LINES = [
    "xtset id t",
    "quietly xtreg x z, re",
    f"scalar drift_u = abs(e(sigma_u) - {SIGMA_U!r})",
    f"scalar drift_e = abs(e(sigma_e) - {SIGMA_E!r})",
    f"scalar su = {SIGMA_U!r}",
    f"scalar se = {SIGMA_E!r}",
    "bysort id: gen double Ti = _N",
    "gen double theta = 1 - se / sqrt(se^2 + Ti * su^2)",
    "quietly summarize theta",
    "scalar drift_min = abs(r(min) - e(thta_min))",
    "scalar drift_max = abs(r(max) - e(thta_max))",
    "bysort id: egen double mx = mean(x)",
    "bysort id: egen double mz = mean(z)",
    "gen double xq = x - theta * mx",
    "gen double zq = z - theta * mz",
    "sort id t",
    export("transforms_quasi", "data", 'st_data(., ("id", "xq", "zq", "theta"))'),
    export(
        "transforms_quasi",
        "sigma_drift",
        '(st_numscalar("drift_u") \\ st_numscalar("drift_e") \\ '
        'st_numscalar("drift_min") \\ st_numscalar("drift_max"))',
    ),
]

case("transforms_quasi", domain=DOMAIN, dataset="sim_re_effects", stata=_QUASI_LINES, atol=1e-9)(
    _quasi
)

# ---------------------------------------------------------------- two-way (D.4.3)

TWOWAY_PANELS = ("balanced", "gaps", "entry_exit", "weak", "delta2", "quarterly", "re_effects")


def _twoway(name: str) -> dict[str, np.ndarray]:
    frame, panel = _declared(name)
    out = panel.within(["x", "z"], effects="twoway")
    rows = panel.index.source_row
    table = np.column_stack(
        [frame["id"].to_numpy(float)[rows], out["x"].to_numpy()[rows], out["z"].to_numpy()[rows]]
    )
    return {"data": table}


def _register_twoway(name: str) -> None:
    tag = f"transforms_twoway_{name}"
    lines = [
        f"xtset id t{_option(name)}",
        "quietly regress x i.id i.t",
        "predict double ex, residuals",
        "quietly regress z i.id i.t",
        "predict double ez, residuals",
        "sort id t",
        export(tag, "data", 'st_data(., ("id", "ex", "ez"))'),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines, atol=1e-8, rtol=1e-8)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        return _twoway(name)


for _name in TWOWAY_PANELS:
    _register_twoway(_name)

# ---------------------------------------------------------------- forward orthogonal deviations (D.4.4)

FOD_PANELS = (
    "balanced",
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "quarterly",
    "unsorted_ids",
    "one_obs_unit",
)

#: The forward deviations written in Stata with running sums taken from the end of each unit.
_FOD_STATA = [
    "gsort id -t",
    "foreach v in x z {",
    "    quietly generate byte ok_`v' = !missing(`v')",
    "    quietly generate double v0_`v' = cond(ok_`v', `v', 0)",
    "    quietly by id: generate double cs_`v' = sum(v0_`v')",
    "    quietly by id: generate double cn_`v' = sum(ok_`v')",
    "    quietly generate double ln_`v' = cn_`v' - ok_`v'",
    "    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * "
    "(`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0",
    "}",
    "sort id t",
]


def _fod(name: str) -> dict[str, np.ndarray]:
    frame, panel = _declared(name)
    out = panel.fod(["x", "z"])
    rows = panel.index.source_row
    table = np.column_stack(
        [frame["id"].to_numpy(float)[rows], out["x"].to_numpy()[rows], out["z"].to_numpy()[rows]]
    )
    return {"data": table}


def _register_fod(name: str) -> None:
    tag = f"transforms_fod_{name}"
    lines = [
        f"xtset id t{_option(name)}",
        *_FOD_STATA,
        export(tag, "data", 'st_data(., ("id", "fod_x", "fod_z"))'),
    ]

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines)
    def _run(name: str = name) -> dict[str, np.ndarray]:
        return _fod(name)


for _name in FOD_PANELS:
    _register_fod(_name)
