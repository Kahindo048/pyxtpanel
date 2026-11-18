"""Hadri test against ``xtunitroot hadri`` in Stata 17 (domain ``diagnostics_unitroot``).

The cases cover the constant and the panel trends, the removal of cross-sectional means, the
heteroskedasticity-robust statistic and the three long-run variance kernels (Bartlett, Parzen and
quadratic-spectral, the last with a fractional bandwidth), on panels of random walks, stationary
series and series with a trend.  Stata returns z, its p-value and the mean and variance of the LM
statistic under the null; all are compared with the number of panels and of periods.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.unitroot.hadri import hadri
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_unitroot"

#: tag -> (dataset, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, str, dict[str, Any]]] = {
    "ha_constant": ("ur_walk", "", {}),
    "ha_trend": ("ur_trend", "trend", {"deterministics": "trend"}),
    "ha_demean": ("ur_mixed", "demean", {"demean": True}),
    "ha_robust": ("ur_stationary", "robust", {"robust": True}),
    "ha_robust_trend": (
        "ur_trend",
        "trend robust demean",
        {"deterministics": "trend", "robust": True, "demean": True},
    ),
    "ha_bartlett": ("ur_mixed", "kernel(bartlett 2)", {"kernel": "bartlett", "lags": 2}),
    "ha_parzen_trend": (
        "ur_trend",
        "trend kernel(parzen 3)",
        {"deterministics": "trend", "kernel": "parzen", "lags": 3},
    ),
    "ha_qs": (
        "ur_long",
        "kernel(quadraticspectral 4.5)",
        {"kernel": "quadraticspectral", "lags": 4.5},
    ),
    "ha_kernel_default": (
        "ur_walk",
        "kernel(bartlett) demean",
        {"kernel": "bartlett", "demean": True},
    ),
}


def _stata(tag: str, options: str) -> list[str]:
    after = f", {options}" if options else ""
    return [
        "quietly xtset id t",
        f"quietly xtunitroot hadri y{after}",
        "scalar hz = r(z)",
        "scalar hp = r(p)",
        "scalar hmu = r(mu)",
        "scalar hvar = r(var)",
        "scalar hng = r(N_g)",
        "scalar hnt = r(N_t)",
        export_scalars(tag, "h", ["hz", "hp", "hmu", "hvar", "hng", "hnt"]),
    ]


def _register(tag: str, dataset: str, options: str, kwargs: dict[str, Any]) -> None:
    def run() -> dict[str, np.ndarray]:
        result = hadri(xtset(panels.load(dataset), "id", "t"), "y", **kwargs)
        assert result.details is not None
        _, mean, variance = result.details["value"]
        return {
            "h": np.array(
                [
                    result.statistic,
                    result.pvalue,
                    mean,
                    variance,
                    result.n_groups,
                    result.n_obs // result.n_groups,
                ],
                dtype=float,
            )
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, options),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _options, _kwargs) in SPECS.items():
    _register(_tag, _dataset, _options, _kwargs)
