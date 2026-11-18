"""Westerlund test against ``xtcointtest westerlund`` in Stata 17 (domain
``diagnostics_cointegration``).

Stata returns the statistic ``r(stat)`` (the standardised variance ratio), its p-value (printed as
``1 - normal(|z|)``, so the left tail ``normal(z)`` is computed in Stata for the comparison), the
number of panels and of observations.  The cases cover the two alternatives (``somepanels``,
``allpanels``), a constant and a trend, one to seven regressors, the removal of cross-sectional
means, a cointegrated panel and an unbalanced panel with very short units.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.cointegration.westerlund import westerlund
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_cointegration"

X7 = [f"x{j}" for j in range(1, 8)]

#: tag -> (dataset, regressors, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, list[str], str, dict[str, Any]]] = {
    "we_some_m1": ("coint_no", ["x"], "", {}),
    "we_all_m1": ("coint_no", ["x"], "allpanels", {"alternative": "all"}),
    "we_some_m3_trend": (
        "coint_wide",
        X7[:3],
        "trend",
        {"deterministics": "trend"},
    ),
    "we_all_m3_trend": (
        "coint_wide",
        X7[:3],
        "trend allpanels",
        {"deterministics": "trend", "alternative": "all"},
    ),
    "we_some_m7": ("coint_wide", X7, "", {}),
    "we_all_m7_trend": (
        "coint_wide",
        X7,
        "trend allpanels",
        {"deterministics": "trend", "alternative": "all"},
    ),
    "we_demean": ("coint_yes", ["x", "z"], "demean", {"demean": True}),
    "we_coint": ("coint_yes", ["x", "z"], "", {}),
    "we_unbalanced": ("coint_unbalanced", ["x", "z"], "", {}),
    "we_unbalanced_all_trend": (
        "coint_unbalanced",
        ["x", "z"],
        "trend allpanels",
        {"deterministics": "trend", "alternative": "all"},
    ),
    "we_unbalanced_demean_all": (
        "coint_unbalanced",
        ["x"],
        "demean allpanels",
        {"demean": True, "alternative": "all"},
    ),
}


def _stata(tag: str, regressors: list[str], options: str) -> list[str]:
    command = f"quietly xtcointtest westerlund y {' '.join(regressors)}"
    if options:
        command += f", {options}"
    # the ado stores the residuals in a variable of the default type (float): double makes the
    # computation the one of the formulas, as for the Kao cases
    return [
        "quietly xtset id t",
        "set type double",
        command,
        "scalar st = r(stat)",
        "scalar pl = normal(r(stat))",
        "scalar ng = r(N_g)",
        "scalar nn = r(N)",
        "set type float",
        export_scalars(tag, "w", ["st", "pl", "ng", "nn"]),
    ]


def _register(
    tag: str, dataset: str, regressors: list[str], options: str, kwargs: dict[str, Any]
) -> None:
    def run() -> dict[str, np.ndarray]:
        result = westerlund(xtset(panels.load(dataset), "id", "t"), "y", regressors, **kwargs)
        assert result.n_groups is not None
        return {
            "w": np.array(
                [result.statistic, result.pvalue, result.n_groups, result.n_obs], dtype=float
            )
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, regressors, options),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _regressors, _options, _kwargs) in SPECS.items():
    _register(_tag, _dataset, _regressors, _options, _kwargs)
