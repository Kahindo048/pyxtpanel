"""Hausman against ``hausman`` in Stata 17 (domain ``diagnostics_postestimation``).

Stata fits ``xtreg, fe`` and ``xtreg, re``, runs ``hausman`` and exports ``e(b)`` and ``e(V)`` of
both fits; the Python side gives those numbers to the function (the constant is left out, which
is what ``hausman`` does by default) and compares ``r(chi2)``, ``r(df)`` and ``r(p)``.
"""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.postestimation.hausman import hausman
from validation.cases._frozen import stata_objects
from validation.registry import case, export, export_scalars

DOMAIN = "diagnostics_postestimation"


def _stata(tag: str, regressors: str, command: str, *, with_p: bool = True) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly xtreg y {regressors}, fe",
        "estimates store fe",
        export(tag, "b_fe", 'st_matrix("e(b)")'),
        export(tag, "v_fe", 'st_matrix("e(V)")'),
        f"quietly xtreg y {regressors}, re",
        "estimates store re",
        export(tag, "b_re", 'st_matrix("e(b)")'),
        export(tag, "v_re", 'st_matrix("e(V)")'),
        f"quietly {command}",
        "scalar hc = r(chi2)",
        "scalar hdf = r(df)",
        "scalar hp = r(p)",
        export_scalars(tag, "h", ["hc", "hdf", "hp"] if with_p else ["hc", "hdf"]),
    ]


def _python(tag: str, *, swapped: bool = False, df: int | None = None) -> dict[str, np.ndarray]:
    s = stata_objects(DOMAIN, tag)
    first, second = ("re", "fe") if swapped else ("fe", "re")
    k = s["b_fe"].size - 1  # leave the constant out, as hausman does
    result = hausman(
        s[f"b_{first}"].ravel(),
        s[f"b_{second}"].ravel(),
        s[f"v_{first}"],
        s[f"v_{second}"],
        keep=range(k),
        df=df,
    )
    shown = [result.statistic, result.df[0]] + ([] if swapped else [result.pvalue])
    out = {"h": np.array(shown)}
    out.update({name: s[name] for name in ("b_fe", "v_fe", "b_re", "v_re")})
    return out


def _register(
    tag: str,
    dataset: str,
    regressors: str,
    command: str,
    *,
    swapped: bool = False,
    df: int | None = None,
) -> None:
    def run() -> dict[str, np.ndarray]:
        return _python(tag, swapped=swapped, df=df)

    case(
        tag,
        domain=DOMAIN,
        dataset=dataset,
        stata=_stata(tag, regressors, command, with_p=not swapped),
        atol=1e-9,
        rtol=1e-9,
    )(run)


_register("hausman_fe_re", "sim_iv_valid", "x z1 z2", "hausman fe re")
_register("hausman_fe_re_invalid", "sim_iv_invalid", "x z1 z2 z3", "hausman fe re")
# Roles swapped: V(re) - V(fe) is negative semi-definite, so Stata reports chi2 < 0 and prints no
# p-value (its r(p) is the upper tail at a negative value, 1, which is not a p-value): the case
# compares the statistic and the degrees of freedom, and the package returns NaN as the p-value.
_register("hausman_swapped_negative", "sim_iv_valid", "x z1 z2", "hausman re fe", swapped=True)
_register("hausman_df_option", "sim_iv_valid", "x z1 z2", "hausman fe re, df(2)", df=2)
