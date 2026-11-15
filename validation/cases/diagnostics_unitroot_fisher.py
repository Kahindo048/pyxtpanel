"""Fisher-type test against ``xtunitroot fisher, dfuller`` in Stata 17 (domain ``diagnostics_unitroot``).

For each case Stata returns the four statistics P, Z, L and Pm with their degrees of freedom and
p-values, the number of panels and of observations, and (through the undocumented ``generate()``
option) the p-value of the Dickey-Fuller test of every unit, floored at the machine epsilon as in
Stata.  All of them are compared with the four runs of ``fisher`` and its per-unit table.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.unitroot.fisher import fisher
from validation import panels
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_unitroot"

#: tag -> (dataset, variable, lags, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, str, int, str, dict[str, Any]]] = {
    "fi_lags1": ("ur_mixed", "y", 1, "", {}),
    "fi_trend_demean": (
        "ur_mixed",
        "y",
        1,
        "trend demean",
        {"deterministics": "trend", "demean": True},
    ),
    "fi_drift": ("ur_mixed", "y", 2, "drift", {"deterministics": "drift"}),
    "fi_entry_exit": ("entry_exit", "x", 0, "demean", {"demean": True}),
    "fi_gaps": ("gaps", "x", 0, "", {}),
    "fi_short_units": ("one_obs_unit", "x", 0, "", {}),
}
METHODS = ("inverse_chi2", "inverse_normal", "logit", "modified")
EPSILON = 2.220446049250313e-16


def _stata(tag: str, variable: str, lags: int, options: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly xtunitroot fisher {variable}, dfuller lags({lags}) {options} generate(pu)",
        "scalar fp = r(P)",
        "scalar fdfp = r(df_P)",
        "scalar fpp = r(p_P)",
        "scalar fz = r(Z)",
        "scalar fpz = r(p_Z)",
        "scalar fl = r(L)",
        "scalar fdfl = r(df_L)",
        "scalar fpl = r(p_L)",
        "scalar fpm = r(Pm)",
        "scalar fppm = r(p_Pm)",
        "scalar fng = r(N_g)",
        "scalar fn = r(N)",
        "preserve",
        "quietly keep if !missing(pu)",
        export(tag, "unit_p", 'st_data(., "pu")'),
        "restore",
        export_scalars(
            tag,
            "f",
            ["fp", "fdfp", "fpp", "fz", "fpz", "fl", "fdfl", "fpl", "fpm", "fppm", "fng", "fn"],
        ),
    ]


def _register(
    tag: str, dataset: str, variable: str, lags: int, options: str, kwargs: dict[str, Any]
) -> None:
    def run() -> dict[str, np.ndarray]:
        panel = xtset(panels.load(dataset), "id", "t")
        r = {m: fisher(panel, variable, lags=lags, method=m, **kwargs) for m in METHODS}  # type: ignore[arg-type]
        details = r["inverse_chi2"].details
        assert details is not None
        p, z, ell, pm = r["inverse_chi2"], r["inverse_normal"], r["logit"], r["modified"]
        return {
            "f": np.array(
                [
                    p.statistic,
                    p.df[0],
                    p.pvalue,
                    z.statistic,
                    z.pvalue,
                    ell.statistic,
                    ell.df[0],
                    ell.pvalue,
                    pm.statistic,
                    pm.pvalue,
                    p.n_groups,
                    p.n_obs,
                ]
            ),
            "unit_p": np.maximum(details["pvalue"].to_numpy(), EPSILON).reshape(-1, 1),
        }

    case(
        tag,
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=_stata(tag, variable, lags, options),
        atol=1e-9,
        rtol=1e-9,
    )(run)


for _tag, (_dataset, _variable, _lags, _options, _kwargs) in SPECS.items():
    _register(_tag, _dataset, _variable, _lags, _options, _kwargs)
