"""Im-Pesaran-Shin test against ``xtunitroot ips`` in Stata 17 (domain ``diagnostics_unitroot``).

Without lags Stata returns t-bar, t-tilde-bar, Z-t-bar and Z-t-tilde-bar with its p-value (missing
when a unit is too short) and, for a balanced panel, the exact critical values; with lags it
returns W-t-bar, its p-value and the average number of lags.  The counts are compared too.  The
short balanced panel has 8 periods (first rows of the moment tables), the unbalanced ones 13 to 23
periods (interpolated columns) and 4 to 10 periods (no p-value), the long one 60 periods.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.unitroot.im_pesaran_shin import im_pesaran_shin
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_unitroot"

#: tag -> (dataset, variable, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, str, str, dict[str, Any]]] = {
    "ips_balanced": ("ur_mixed", "y", "", {}),
    "ips_short": ("balanced", "x", "", {}),
    "ips_trend": ("ur_trend", "y", "trend", {"deterministics": "trend"}),
    "ips_demean": ("ur_walk", "y", "demean", {"demean": True}),
    "ips_long": ("ur_long", "y", "", {}),
    "ips_unbalanced": ("ur_unbalanced", "y", "", {}),
    "ips_unbalanced_short": ("entry_exit", "x", "", {}),
    "ips_lags1": ("ur_mixed", "y", "lags(1)", {"lags": 1}),
    "ips_lags_aic": ("ur_mixed", "y", "lags(aic 3)", {"lags": 3, "lag_selection": "aic"}),
    "ips_lags_trend_unbalanced": (
        "ur_unbalanced",
        "y",
        "lags(2) trend",
        {"lags": 2, "deterministics": "trend"},
    ),
    "ips_lags_hqic_long": (
        "ur_long",
        "y",
        "lags(hqic 4) trend",
        {"lags": 4, "lag_selection": "hqic", "deterministics": "trend"},
    ),
    "ips_nointerpolate": (
        "ur_unbalanced",
        "y",
        "lags(1) nointerpolate",
        {"lags": 1, "interpolate": False},
    ),
}

NO_LAGS = ("zttildebar", "p_zttildebar", "tbar", "ttildebar", "zt", "cv_1", "cv_5", "cv_10")


def _stata(tag: str, variable: str, options: str) -> list[str]:
    after = f", {options}" if options else ""
    lines = ["quietly xtset id t", f"quietly xtunitroot ips {variable}{after}"]
    names = []
    if "lags" not in options:
        for key in NO_LAGS:
            lines.append(f"scalar i_{key} = r({key})")
            names.append(f"i_{key}")
    else:
        selected = "aic" in options or "hqic" in options
        lines += ["scalar i_wtbar = r(wtbar)", "scalar i_p_wtbar = r(p_wtbar)"]
        lines.append(f"scalar i_lags = r({'lagm' if selected else 'lags'})")
        names += ["i_wtbar", "i_p_wtbar", "i_lags"]
    for key in ("N", "N_g", "N_t"):
        lines.append(f"scalar i_{key} = r({key})")
        names.append(f"i_{key}")
    lines.append(export_scalars(tag, "ips", names))
    return lines


def _register(tag: str, dataset: str, variable: str, options: str, kwargs: dict[str, Any]) -> None:
    def run() -> dict[str, np.ndarray]:
        result = im_pesaran_shin(xtset(panels.load(dataset), "id", "t"), variable, **kwargs)
        assert result.details is not None
        d = result.details["value"]
        if "lags" not in kwargs:
            head = [
                result.statistic,
                result.pvalue,
                d["t_bar"],
                d["t_tilde_bar"],
                d["z_t_bar"],
                d["cv_1"],
                d["cv_5"],
                d["cv_10"],
            ]
        else:
            head = [result.statistic, result.pvalue, d["lags"]]
        counts = [d["n_equations"], result.n_groups, d["average_periods"]]
        return {"ips": np.array([*head, *counts])}

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, variable, options),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _variable, _options, _kwargs) in SPECS.items():
    _register(_tag, _dataset, _variable, _options, _kwargs)
