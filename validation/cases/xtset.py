"""``xtset`` against Stata (C.6.1): time range, delta, gaps, balance and counts.

For every simulated panel Stata runs ``xtset`` and the do-file exports ``r(tmin)``, ``r(tmax)``,
``r(tdelta)``, ``r(gaps)``, ``r(imin)``, ``r(imax)``, the balance class (1 strongly, 2 weakly,
3 unbalanced) and the counts of observations, units and distinct periods.  The refused panels
check the return code.
"""

from __future__ import annotations

import numpy as np

from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "xtset"
NAMES = ["tmin", "tmax", "tdelta", "gaps", "imin", "imax", "bal", "nobs", "nunits", "nperiods"]

_COUNTS = [
    "scalar nobs = _N",
    "preserve",
    "drop if id >= . | t >= .",
    "scalar nobs = _N",
    "bysort id: gen byte _first = (_n == 1)",
    "quietly summarize _first",
    "scalar nunits = r(sum)",
    "bysort t: gen byte _firstt = (_n == 1)",
    "quietly summarize _firstt",
    "scalar nperiods = r(sum)",
    "restore",
]


def _stata_lines(tag: str, delta: int) -> list[str]:
    option = f", delta({delta})" if delta != 1 else ""
    return [
        f"xtset id t{option}",
        "scalar tmin = r(tmin)",
        "scalar tmax = r(tmax)",
        "scalar tdelta = r(tdelta)",
        "scalar gaps = r(gaps)",
        "scalar imin = r(imin)",
        "scalar imax = r(imax)",
        'scalar bal = cond("`r(balanced)\'" == "strongly balanced", 1, '
        'cond("`r(balanced)\'" == "weakly balanced", 2, 3))',
        *_COUNTS,
        export_scalars(tag, "xtset", NAMES),
    ]


def _python(name: str, missing_keys: str = "error") -> dict[str, np.ndarray]:
    from pyxtpanel import xtset

    spec = panels.SPECS[name]
    panel = xtset(
        panels.load(name),
        "id",
        "t",
        delta=spec.delta,
        missing_keys=missing_keys,  # type: ignore[arg-type]
    )
    rank = panel.index.rank
    first, last = int(rank.min()), int(rank.max())
    labels = panel.unit_labels.astype(float)
    code = {"strongly balanced": 1, "weakly balanced": 2, "unbalanced": 3}[panel.balanced]
    values = [
        float(panel.axis.stata_values(first)[0]),
        float(panel.axis.stata_values(last)[0]),
        float(panel.delta),
        float(panel.has_gaps),
        float(labels.min()),
        float(labels.max()),
        float(code),
        float(panel.n_obs),
        float(panel.n_units),
        float(panel.n_periods),
    ]
    return {"xtset": np.array(values).reshape(-1, 1)}


def _register_panel(name: str, missing_keys: str = "error") -> None:
    spec = panels.SPECS[name]
    tag = f"xtset_{name}"

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=_stata_lines(tag, spec.delta))
    def _run(name: str = name, missing_keys: str = missing_keys) -> dict[str, np.ndarray]:
        return _python(name, missing_keys)


for _name in (
    "balanced",
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "quarterly",
    "monthly",
    "unsorted_ids",
    "constants",
    "one_obs_unit",
    "single_unit",
    "single_period",
    "all_missing_var",
):
    _register_panel(_name)
_register_panel("missing_keys", missing_keys="drop")


def _refused(name: str) -> None:
    tag = f"xtset_refused_{name}"

    @case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{name}",
        stata=["capture xtset id t", "scalar rc = _rc", export_scalars(tag, "rc", ["rc"])],
    )
    def _run(name: str = name) -> dict[str, np.ndarray]:
        from pyxtpanel import PyxtpanelError, xtset

        try:
            xtset(panels.HOSTILE_REFUSED[name], "id", "t")
        except PyxtpanelError as error:
            return {"rc": np.array([[float(error.stata_code or 0)]])}
        return {"rc": np.array([[0.0]])}


for _name in panels.HOSTILE_REFUSED:
    _refused(_name)
