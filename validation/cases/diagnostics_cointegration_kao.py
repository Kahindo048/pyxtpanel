"""Kao tests against ``xtcointtest kao`` in Stata 17 (domain ``diagnostics_cointegration``).

Each configuration is registered twice, because ``xtcointtest`` hands its ``t`` ratios to Mata
through ``strofreal``, which keeps seven significant digits: the statistics built on them (DF t,
ADF t and the unadjusted DF t) agree with this package to 1e-6, the others (the modified DF t, the
unadjusted modified DF t, the lags and the counts) to 1e-9.  ``set type double`` removes the
single-precision storage of the residuals that ``xtcointtest`` creates.  The p-values are
``normal(statistic)`` computed in Stata: ``xtcointtest`` reports ``1 - normal(abs(statistic))``,
which differs for a positive statistic.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.cointegration.kao import kao
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_cointegration"

#: tag -> (dataset, regressors, Stata options, Python keyword arguments)
SPECS: dict[str, tuple[str, tuple[str, ...], str, dict[str, Any]]] = {
    "kao_default": ("coint_yes", ("x", "z"), "", {}),
    "kao_nocoint": ("coint_no", ("x", "z"), "", {}),
    "kao_bartlett3": ("coint_yes", ("x", "z"), "kernel(bartlett 3)", {"kernel_lags": 3}),
    "kao_parzen": (
        "coint_yes",
        ("x", "z"),
        "kernel(parzen 2) lags(2)",
        {"kernel": "parzen", "kernel_lags": 2, "adf_lags": 2},
    ),
    "kao_parzen_auto": ("coint_no", ("x", "z"), "kernel(parzen nwest)", {"kernel": "parzen"}),
    "kao_adf2": ("coint_yes", ("x", "z"), "lags(2)", {"adf_lags": 2}),
    "kao_adf_aic": (
        "coint_no",
        ("x", "z"),
        "lags(aic 3)",
        {"adf_lags": 3, "adf_lag_selection": "aic"},
    ),
    "kao_adf_bic": (
        "coint_yes",
        ("x", "z"),
        "lags(bic 3)",
        {"adf_lags": 3, "adf_lag_selection": "bic"},
    ),
    "kao_adf_hqic": (
        "coint_no",
        ("x", "z"),
        "lags(hqic 4)",
        {"adf_lags": 4, "adf_lag_selection": "hqic"},
    ),
    "kao_demean": ("coint_yes", ("x", "z"), "demean", {"demean": True}),
    "kao_one_regressor": ("coint_no", ("x",), "", {}),
    "kao_qs": (
        "coint_yes",
        ("x", "z"),
        "kernel(quadraticspectral 3)",
        {"kernel": "quadraticspectral", "kernel_lags": 3},
    ),
    "kao_qs_auto": (
        "coint_no",
        ("x", "z"),
        "kernel(quadraticspectral nwest)",
        {"kernel": "quadraticspectral"},
    ),
}
#: Stata column order of ``r(stats)`` and the names used by the package
ORDER = ("modified_df", "df", "adf", "unadjusted_modified_df", "unadjusted_df")
EXACT = ("modified_df", "unadjusted_modified_df")
ROUNDED = ("df", "adf", "unadjusted_df")


def _stata(variables: str, options: str) -> list[str]:
    after = f", {options}" if options else ""
    lines = [
        "quietly xtset id t",
        "set type double",
        f"quietly xtcointtest kao y {variables}{after}",
        "matrix S = r(stats)",
    ]
    for i, name in enumerate(ORDER, start=1):
        lines += [f"scalar s_{name} = S[1, {i}]", f"scalar p_{name} = normal(S[1, {i}])"]
    lines += [
        "scalar k_lags = r(hac_lagm)",
        "scalar k_adf = r(adf_lags)",
        "scalar k_ng = r(N_g)",
        "scalar k_nt = r(N_t)",
        "scalar k_n = r(N)",
        "set type float",
    ]
    return lines


def _register(
    tag: str, dataset: str, variables: tuple[str, ...], options: str, kwargs: dict[str, Any]
) -> None:
    def run() -> dict[str, np.ndarray]:
        result = kao(xtset(panels.load(dataset), "id", "t"), "y", list(variables), **kwargs)
        assert result.details is not None
        table = result.details
        return {
            "exact": np.array(
                [table.loc[k, "statistic"] for k in EXACT]
                + [table.loc[k, "pvalue"] for k in EXACT]
                + [
                    float(result.options["kernel_lags"]),  # type: ignore[arg-type]
                    float(result.options["adf_lags"]),  # type: ignore[arg-type]
                    result.n_groups,
                    result.n_obs // result.n_groups,
                    result.n_obs,
                ]
            ),
            "rounded": np.array(
                [table.loc[k, "statistic"] for k in ROUNDED]
                + [table.loc[k, "pvalue"] for k in ROUNDED]
            ),
        }

    def stata_exact() -> list[str]:
        names = [f"s_{k}" for k in EXACT] + [f"p_{k}" for k in EXACT]
        names += ["k_lags", "k_adf", "k_ng", "k_nt", "k_n"]
        return [*_stata(" ".join(variables), options), export_scalars(tag, "exact", names)]

    def stata_rounded() -> list[str]:
        names = [f"s_{k}" for k in ROUNDED] + [f"p_{k}" for k in ROUNDED]
        return [*_stata(" ".join(variables), options), export_scalars(f"{tag}_t", "rounded", names)]

    def run_exact() -> dict[str, np.ndarray]:
        return {"exact": run()["exact"]}

    def run_rounded() -> dict[str, np.ndarray]:
        return {"rounded": run()["rounded"]}

    case(tag, domain=DOMAIN, dataset=f"sim_{dataset}", stata=stata_exact(), atol=1e-9, rtol=1e-9)(
        run_exact
    )
    case(
        f"{tag}_t",
        domain=DOMAIN,
        dataset=f"sim_{dataset}",
        stata=stata_rounded(),
        atol=1e-6,
        rtol=1e-6,
    )(run_rounded)


for _tag, (_dataset, _variables, _options, _kwargs) in SPECS.items():
    _register(_tag, _dataset, _variables, _options, _kwargs)
