"""Pedroni test against ``xtcointtest pedroni`` in Stata 17 (domain ``diagnostics_cointegration``).

Stata returns the standardised statistics in ``r(stats)`` (the group statistics with the default
panel-specific autoregressive parameter, the four panel statistics with ``ar(same)``), the number
of augmented lags, of panels and of observations.  Its p-values are printed as ``1 - normal(|z|)``;
the one-sided tail of each statistic (left for the six others, right for the variance ratio) is
computed in Stata for the comparison.  The cases fix the number of lags of the kernel (Stata's
automatic bandwidth is not reproduced) and cover the three deterministic specifications, one to
seven regressors, the three kernels, the augmented lags and their choice by information criteria,
the removal of cross-sectional means, a cointegrated panel and unbalanced panels.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pyxtpanel import xtset
from pyxtpanel.diagnostics.cointegration.pedroni import pedroni
from validation import panels
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_cointegration"

X7 = [f"x{j}" for j in range(1, 8)]
SAME = {"ar": "same"}

#: tag -> (dataset, regressors, Stata options, Python keyword arguments); the number of lags of the
#: kernel is 2 unless the options say otherwise
SPECS: dict[str, tuple[str, list[str], str, dict[str, Any]]] = {
    "pe_group_m1": ("coint_no", ["x"], "", {}),
    "pe_panel_m1": ("coint_no", ["x"], "ar(same)", SAME),
    "pe_group_m3_trend": ("coint_wide", X7[:3], "trend", {"deterministics": "trend"}),
    "pe_panel_m3_trend": (
        "coint_wide",
        X7[:3],
        "trend ar(same)",
        {"deterministics": "trend", **SAME},
    ),
    "pe_group_m2_none": (
        "coint_yes",
        ["x", "z"],
        "noconstant",
        {"deterministics": "none"},
    ),
    "pe_panel_m2_none": (
        "coint_yes",
        ["x", "z"],
        "noconstant ar(same)",
        {"deterministics": "none", **SAME},
    ),
    "pe_group_m7": ("coint_wide", X7, "", {}),
    "pe_panel_m7_trend": (
        "coint_wide",
        X7,
        "trend ar(same)",
        {"deterministics": "trend", **SAME},
    ),
    "pe_parzen": (
        "coint_yes",
        ["x", "z"],
        "kernel(parzen 3) ar(same)",
        {"kernel": "parzen", "kernel_lags": 3, **SAME},
    ),
    "pe_quadraticspectral": (
        "coint_yes",
        ["x", "z"],
        "kernel(quadraticspectral 3)",
        {"kernel": "quadraticspectral", "kernel_lags": 3},
    ),
    "pe_kernel0": (
        "coint_no",
        ["x", "z"],
        "kernel(bartlett 0) ar(same)",
        {"kernel_lags": 0, **SAME},
    ),
    "pe_adf2": ("coint_yes", ["x", "z"], "lags(2)", {"adf_lags": 2}),
    "pe_adf_aic": (
        "coint_yes",
        ["x", "z"],
        "lags(aic 3)",
        {"adf_lags": 3, "adf_lag_selection": "aic"},
    ),
    "pe_adf_bic": (
        "coint_wide",
        X7[:2],
        "lags(bic 3) ar(same)",
        {"adf_lags": 3, "adf_lag_selection": "bic", **SAME},
    ),
    "pe_adf_hqic": (
        "coint_no",
        ["x", "z"],
        "lags(hqic 4)",
        {"adf_lags": 4, "adf_lag_selection": "hqic"},
    ),
    "pe_demean": ("coint_yes", ["x", "z"], "demean", {"demean": True}),
    "pe_unbalanced": ("coint_unbalanced", ["x", "z"], "ar(same)", SAME),
    "pe_unbalanced_panel_trend": (
        "coint_unbalanced",
        ["x", "z"],
        "trend",
        {"deterministics": "trend"},
    ),
    "pe_unbalanced_none": (
        "coint_unbalanced",
        ["x"],
        "noconstant ar(same)",
        {"deterministics": "none", **SAME},
    ),
    "pe_coint": ("coint_yes", ["x", "z"], "ar(same)", SAME),
}


def _stata(tag: str, regressors: list[str], options: str) -> list[str]:
    same = "ar(same)" in options
    kernel = "" if "kernel(" in options else "kernel(bartlett 2)"
    command = f"quietly xtcointtest pedroni y {' '.join(regressors)}, {kernel} {options}"
    names = ["nu", "prho", "ptpp", "padf"] if same else ["grho", "gtpp", "gadf"]
    lines = [
        "quietly xtset id t",
        "set type double",
        command.replace(",  ", ", ").replace(", ,", ","),
        "matrix S = r(stats)",
    ]
    scalars = []
    for index, name in enumerate(names, start=1):
        lines.append(f"scalar z_{name} = S[1,{index}]")
        tail = f"1 - normal(z_{name})" if name == "nu" else f"normal(z_{name})"
        lines.append(f"scalar p_{name} = {tail}")
        scalars += [f"z_{name}", f"p_{name}"]
    lines += [
        "scalar nl = r(adf_lags)",
        "scalar ng = r(N_g)",
        "scalar nn = r(N)",
        "set type float",
        export_scalars(tag, "p", [*scalars, "nl", "ng", "nn"]),
    ]
    return lines


def _register(
    tag: str, dataset: str, regressors: list[str], options: str, kwargs: dict[str, Any]
) -> None:
    same = kwargs.get("ar") == "same"
    keys = ("nu", "rho", "pp", "adf") if same else ("rho", "pp", "adf")

    def run() -> dict[str, np.ndarray]:
        options_py = {"kernel_lags": 2, **kwargs}
        panel = xtset(panels.load(dataset), "id", "t")
        results = [pedroni(panel, "y", regressors, statistic=key, **options_py) for key in keys]  # type: ignore[arg-type]
        first = results[0]
        assert first.details is not None and first.n_groups is not None
        out: list[float] = []
        for result in results:
            out += [result.statistic, result.pvalue]
        out += [float(first.options["adf_lags"]), float(first.n_groups), float(first.n_obs)]  # type: ignore[arg-type]
        return {"p": np.array(out, dtype=float)}

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
