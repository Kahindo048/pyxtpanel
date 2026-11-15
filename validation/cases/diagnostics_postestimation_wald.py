"""Wald against ``test`` in Stata 17 (domain ``diagnostics_postestimation``).

Stata fits a model, runs ``test`` and exports ``e(b)``, ``e(V)`` and ``e(df_r)``; the Python side
gives those numbers and the restrictions to the function and compares the statistic (F or
chi-square), the degrees of freedom and the p-value.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.postestimation.wald import wald
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_postestimation"


def _selection(k: int, positions: tuple[int, ...]) -> np.ndarray:
    out = np.zeros((len(positions), k))
    out[np.arange(len(positions)), positions] = 1.0
    return out


def _stata(tag: str, fit: str, test: str, *, f_version: bool) -> list[str]:
    stat = "r(F)" if f_version else "r(chi2)"
    lines = [
        "quietly xtset id t",
        f"quietly {fit}",
        export(tag, "b", 'st_matrix("e(b)")'),
        export(tag, "v", 'st_matrix("e(V)")'),
        f"quietly {test}",
        f"scalar ws = {stat}",
        "scalar wdf = r(df)",
        "scalar wp = r(p)",
    ]
    if f_version:
        lines += ["scalar dfr = r(df_r)", export_scalars(tag, "dfr", ["dfr"])]
    return [*lines, export_scalars(tag, "w", ["ws", "wdf", "wp"])]


def _register(
    tag: str,
    dataset: str,
    fit: str,
    test: str,
    positions: tuple[int, ...],
    values: tuple[float, ...],
    *,
    f_version: bool,
) -> None:
    def run() -> dict[str, np.ndarray]:
        s = stata_objects(DOMAIN, tag)
        b, v = s["b"].ravel(), s["v"]
        df_resid = float(s["dfr"].ravel()[0]) if f_version else None
        result = wald(b, v, _selection(b.size, positions), np.array(values), df_resid=df_resid)
        out = {"w": np.array([result.statistic, result.df[0], result.pvalue]), "b": s["b"], "v": v}
        if f_version:
            out["dfr"] = np.array([result.df[1]])
        return out

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, fit, test, f_version=f_version),
        atol=1e-9,
        rtol=1e-9,
    )(run)


# Coefficient order of e(b): the regressors in the order of the command, then the constant.
_register(
    "wald_regress_f",
    "sim_iv_valid",
    "regress y x z1 z2 z3, vce(cluster id)",
    "test x z1",
    (0, 1),
    (0.0, 0.0),
    f_version=True,
)
_register(
    "wald_xtreg_fe_cluster_f",
    "sim_iv_valid",
    "xtreg y x z1 z2, fe vce(cluster id)",
    "test x z2",
    (0, 2),
    (0.0, 0.0),
    f_version=True,
)
_register(
    "wald_ivregress_chi2",
    "sim_iv_valid",
    "ivregress 2sls y (x = z1 z2 z3) z4, vce(cluster id)",
    "test x z4",
    (0, 1),
    (0.0, 0.0),
    f_version=False,
)
_register(
    "wald_nonzero_values",
    "sim_iv_valid",
    "regress y x z1 z2 z3",
    "test (x = 0.8) (z1 = 0.1)",
    (0, 1),
    (0.8, 0.1),
    f_version=True,
)
