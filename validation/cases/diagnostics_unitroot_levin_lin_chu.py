"""Levin-Lin-Chu test against ``xtunitroot llc`` in Stata 17 (domain ``diagnostics_unitroot``).

The long panel (60 periods) exercises the interpolated rows of the adjustment table for the three
models; the short panels (12 and 20 periods) the first row.  Stata returns the unadjusted and
adjusted statistics with every intermediate quantity (``delta``, its standard error, the
adjustments, the average number of equations, ``sbar``), the average numbers of ADF lags and of
kernel lags used, and the counts; all are compared.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.unitroot.levin_lin_chu import levin_lin_chu
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_unitroot"

#: tag -> (dataset, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, str, dict[str, Any]]] = {
    "llc_lags1": ("ur_long", "lags(1)", {"lags": 1}),
    "llc_none": ("ur_long", "noconstant lags(1)", {"deterministics": "none", "lags": 1}),
    "llc_none_aic": (
        "ur_walk",
        "noconstant lags(aic 3)",
        {"deterministics": "none", "lags": 3, "lag_selection": "aic"},
    ),
    "llc_trend": ("ur_long", "trend lags(2)", {"deterministics": "trend", "lags": 2}),
    "llc_demean": ("ur_stationary", "demean lags(1)", {"demean": True, "lags": 1}),
    "llc_aic": ("ur_mixed", "lags(aic 3)", {"lags": 3, "lag_selection": "aic"}),
    "llc_hqic_parzen": (
        "ur_mixed",
        "lags(hqic 3) kernel(parzen 3)",
        {"lags": 3, "lag_selection": "hqic", "kernel": "parzen", "bandwidth": 3},
    ),
    "llc_parzen_nwest": (
        "ur_mixed",
        "lags(2) kernel(parzen nwest)",
        {"lags": 2, "kernel": "parzen", "bandwidth": "nwest"},
    ),
    "llc_bartlett_nwest": (
        "ur_long",
        "lags(1) kernel(bartlett nwest)",
        {"lags": 1, "bandwidth": "nwest"},
    ),
    "llc_qs_fixed": (
        "ur_mixed",
        "lags(1) kernel(andrews 4.5)",
        {"lags": 1, "kernel": "quadraticspectral", "bandwidth": 4.5},
    ),
    "llc_qs_nwest": (
        "ur_mixed",
        "lags(1) kernel(andrews nwest)",
        {"lags": 1, "kernel": "quadraticspectral", "bandwidth": "nwest"},
    ),
    "llc_nointerpolate": (
        "ur_long",
        "lags(1) nointerpolate",
        {"lags": 1, "interpolate": False},
    ),
}

#: Stata return values in the order of the Python vector built below
RETURNS = (
    "tds",
    "td",
    "delta",
    "se_delta",
    "mu_adj",
    "sig_adj",
    "ttilde",
    "sbar",
    "Var_ep",
    "p_tds",
    "N",
    "N_g",
    "N_t",
)


def _stata(tag: str, options: str) -> list[str]:
    lines = ["quietly xtset id t", f"quietly xtunitroot llc y, {options}"]
    names = []
    for key in RETURNS:
        lines.append(f"scalar l_{key} = r({key})")
        names.append(f"l_{key}")
    adf = "adf_lagm" if "aic" in options or "hqic" in options else "adf_lags"
    lines.append(f"scalar l_adf = r({adf})")
    fixed = "parzen 3" in options or "andrews 4.5" in options  # a number of kernel lags was given
    lines.append("scalar l_hac = r(hac_lags)" if fixed else "scalar l_hac = r(hac_lagm)")
    names += ["l_adf", "l_hac"]
    lines.append(export_scalars(tag, "llc", names))
    return lines


def _register(tag: str, dataset: str, options: str, kwargs: dict[str, Any]) -> None:
    def run() -> dict[str, np.ndarray]:
        result = levin_lin_chu(xtset(panels.load(dataset), "id", "t"), "y", **kwargs)
        assert result.details is not None
        d = result.details["value"]
        return {
            "llc": np.array(
                [
                    result.statistic,
                    d["t_delta"],
                    d["delta"],
                    d["se_delta"],
                    d["mu_adj"],
                    d["sigma_adj"],
                    d["t_tilde"],
                    d["s_bar"],
                    d["var_e"],
                    result.pvalue,
                    d["t_tilde"] * result.n_groups,
                    result.n_groups,
                    result.n_obs / result.n_groups,
                    d["adf_lags"],
                    d["hac_lags"],
                ]
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
