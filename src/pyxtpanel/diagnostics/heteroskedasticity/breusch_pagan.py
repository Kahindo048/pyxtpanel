"""Breusch-Pagan / Cook-Weisberg test for heteroskedasticity, from bare arrays (category B, L.8.1).

Validated against Stata 17: ``estat hettest`` after ``regress`` (fitted values, ``rhs``, a
variable list, ``iid`` and ``fstat``) on a homoskedastic and on a heteroskedastic panel; the
residuals and variables are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_heteroskedasticity_breusch_pagan.py``.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, as_array, register_test
from pyxtpanel.diagnostics._shared.ols import ols_fit

_i18n.register(
    "diag",
    {
        "diag.breusch_pagan.title": {
            "en": "Breusch-Pagan / Cook-Weisberg test for heteroskedasticity",
            "fr": "Test d'hétéroscédasticité de Breusch-Pagan / Cook-Weisberg",
        },
        "diag.breusch_pagan.null": {
            "en": "constant variance",
            "fr": "variance constante",
        },
        "diag.breusch_pagan.alternative": {
            "en": "the variance depends on the variables tested",
            "fr": "la variance dépend des variables testées",
        },
        "diag.breusch_pagan.bad_variant": {
            "en": "variant must be 'normal', 'iid' or 'fstat', got {value!r}.",
            "fr": "variant doit valoir 'normal', 'iid' ou 'fstat', reçu {value!r}.",
        },
        "diag.breusch_pagan.rows": {
            "en": "variables must have one row per residual ({n}), got {got}.",
            "fr": "variables doit avoir une ligne par résidu ({n}), reçu {got}.",
        },
        "diag.breusch_pagan.zero_residuals": {
            "en": "all residuals are zero: the variance of the errors cannot be tested.",
            "fr": "tous les résidus sont nuls : la variance des erreurs ne peut pas être testée.",
        },
        "diag.breusch_pagan.no_variation": {
            "en": "the variables tested are constant or collinear with the constant: nothing to "
            "test.",
            "fr": "les variables testées sont constantes ou colinéaires à la constante : rien "
            "à tester.",
        },
        "diag.breusch_pagan.too_few": {
            "en": "{n} observations are too few for {k} independent variables.",
            "fr": "{n} observations sont trop peu pour {k} variables indépendantes.",
        },
    },
)


@register_test(
    "breusch_pagan",
    TestSpec(
        family="heteroskedasticity",
        category="B",
        wave=1,
        reference="Breusch, T. S. and Pagan, A. R. (1979), Econometrica 47(5), 1287-1294; "
        "Cook, R. D. and Weisberg, S. (1983), Biometrika 70(1), 1-10",
        stata="estat hettest",
        validation=(
            "bp_fitted",
            "bp_rhs",
            "bp_variables",
            "bp_iid",
            "bp_fstat",
            "bp_homoskedastic",
        ),
        inputs="arrays",
        balanced="na",
        r="lmtest::bptest",
        summary="Test that the variance of the errors depends on given variables.",
    ),
)
def breusch_pagan(
    residuals: object,
    variables: object,
    *,
    variant: Literal["normal", "iid", "fstat"] = "normal",
) -> TestResult:
    """Return the Breusch-Pagan / Cook-Weisberg test that the variance depends on ``variables``.

    The squared residuals divided by their mean, ``e^2 / (RSS / N)``, are regressed on a constant
    and ``variables``; the statistic is half the explained sum of squares (``"normal"``, the
    Breusch-Pagan form under normal errors, the default of Stata's ``estat hettest``),
    ``N * R^2`` (``"iid"``, Koenker's form that does not need normality) or the regression's
    ``F`` (``"fstat"``).  Pass the fitted values for the default test, the regressors for
    Stata's ``rhs`` or any variables suspected to drive the variance.

    The test is for residuals of a least-squares fit with a constant and homoskedastic-error
    estimation (not robust or clustered, for which Stata refuses it too).

    Parameters
    ----------
    residuals
        ``(n,)`` residuals of the fitted model.
    variables
        ``(n,)`` or ``(n, m)`` variables the variance may depend on, without a constant;
        constant or collinear columns are dropped and not counted in the degrees of freedom.
    variant
        ``"normal"``, ``"iid"`` or ``"fstat"``.

    Raises
    ------
    SpecificationError
        If ``variant`` is unknown or the shapes do not match.
    NumericalError
        If the residuals are all zero or the variables carry no information.
    """
    # Stata: estat hettest -> u = e^2 / (e(rss) / e(N)); regress u on varlist; chi2 = e(mss) / 2
    # (default), e(N) * e(r2) (iid) or e(F) (fstat); df = e(df_m)
    if variant not in ("normal", "iid", "fstat"):
        raise SpecificationError(
            _i18n.message("diag.breusch_pagan.bad_variant", value=variant),
            key="diag.breusch_pagan.bad_variant",
        )
    e = as_array("residuals", residuals, ndim=1)
    raw = np.asarray(variables, dtype=np.float64)
    z = as_array("variables", raw[:, None] if raw.ndim == 1 else raw, ndim=2)
    if z.shape[0] != e.size:
        raise SpecificationError(
            _i18n.message("diag.breusch_pagan.rows", n=e.size, got=z.shape[0]),
            key="diag.breusch_pagan.rows",
        )
    mean_square = float(e @ e) / e.size
    if mean_square == 0.0:
        raise NumericalError(
            _i18n.message("diag.breusch_pagan.zero_residuals"),
            key="diag.breusch_pagan.zero_residuals",
        )
    fit = ols_fit(e**2 / mean_square, z)
    if fit.df_model < 1:
        raise NumericalError(
            _i18n.message("diag.breusch_pagan.no_variation"), key="diag.breusch_pagan.no_variation"
        )
    df: tuple[float, ...]
    if variant == "fstat":
        if fit.df_resid < 1:
            raise NumericalError(
                _i18n.message("diag.breusch_pagan.too_few", n=e.size, k=fit.df_model),
                key="diag.breusch_pagan.too_few",
            )
        statistic, name = fit.f_statistic, "F"
        pvalue = float(stats.f.sf(statistic, fit.df_model, fit.df_resid))
        df = (float(fit.df_model), float(fit.df_resid))
    else:
        statistic = e.size * fit.r2 if variant == "iid" else fit.mss / 2.0
        name, df = "chi2", (float(fit.df_model),)
        pvalue = float(stats.chi2.sf(statistic, fit.df_model))
    return TestResult(
        name="breusch_pagan",
        statistic=statistic if math.isfinite(statistic) else math.nan,
        pvalue=pvalue,
        statistic_name=name,
        df=df,
        n_obs=e.size,
        options={"variant": variant, "variables": fit.df_model},
    )
