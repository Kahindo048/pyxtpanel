"""Breitung test against ``xtunitroot breitung`` in Stata 17 (domain ``diagnostics_unitroot``).

The cases cover the three specifications (constant, no constant, panel trends), prewhitening with
and without lags, the removal of cross-sectional means and the statistic robust to cross-sectional
correlation, which needs more periods than panels (hence the long panel).  Stata returns lambda
and its p-value, and with ``robust`` the robust statistic and its p-value; all are compared with
the number of panels and of periods.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.unitroot.breitung import breitung
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_unitroot"

#: tag -> (dataset, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, str, dict[str, Any]]] = {
    "br_constant": ("ur_walk", "", {}),
    "br_none": ("ur_stationary", "noconstant", {"deterministics": "none"}),
    "br_trend": ("ur_trend", "trend", {"deterministics": "trend"}),
    "br_lags2": ("ur_mixed", "lags(2)", {"lags": 2}),
    "br_trend_lags": ("ur_mixed", "trend lags(1)", {"deterministics": "trend", "lags": 1}),
    "br_none_lags": ("ur_long", "noconstant lags(2)", {"deterministics": "none", "lags": 2}),
    "br_demean": ("ur_mixed", "demean lags(1)", {"demean": True, "lags": 1}),
    "br_robust": ("ur_long", "robust lags(2)", {"robust": True, "lags": 2}),
    "br_robust_trend_lags": (
        "ur_long",
        "trend robust lags(3)",
        {"deterministics": "trend", "robust": True, "lags": 3},
    ),
}


def _stata(tag: str, options: str, robust: bool) -> list[str]:
    scalars = ["blam", "bp", "bng", "bnt"]
    lines = [
        "quietly xtset id t",
        f"quietly xtunitroot breitung y, {options}" if options else "quietly xtunitroot breitung y",
        "scalar blam = r(lambda)",
        "scalar bp = r(p)",
        "scalar bng = r(N_g)",
        "scalar bnt = r(N_t)",
    ]
    if robust:
        lines += ["scalar blr = r(lrobust)", "scalar bpr = r(p_lrobust)"]
        scalars += ["blr", "bpr"]
    return [*lines, export_scalars(tag, "b", scalars)]


def _register(tag: str, dataset: str, options: str, kwargs: dict[str, Any]) -> None:
    def run() -> dict[str, np.ndarray]:
        result = breitung(xtset(panels.load(dataset), "id", "t"), "y", **kwargs)
        assert result.details is not None
        values = result.details["value"]
        out = [
            values["lambda"],
            values["pvalue_lambda"],
            result.n_groups,
            result.n_obs // result.n_groups,
        ]
        if kwargs.get("robust"):
            out += [result.statistic, result.pvalue]
        return {"b": np.array(out, dtype=float)}

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, options, bool(kwargs.get("robust"))),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _options, _kwargs) in SPECS.items():
    _register(_tag, _dataset, _options, _kwargs)
