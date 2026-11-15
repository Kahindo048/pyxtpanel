"""White's test for heteroskedasticity, from bare arrays (category B, L.8.2).

Validated against Stata 17: ``estat imtest, white double`` after ``regress`` on a
homoskedastic and on a heteroskedastic panel, with one to four regressors; the residuals and
regressors are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_heteroskedasticity_white.py``.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, as_array, register_test
from pyxtpanel.diagnostics._shared.ols import ols_fit

_i18n.register(
    "diag",
    {
        "diag.white.title": {
            "en": "White's test for heteroskedasticity",
            "fr": "Test d'hétéroscédasticité de White",
        },
        "diag.white.null": {"en": "homoskedasticity", "fr": "homoscédasticité"},
        "diag.white.alternative": {
            "en": "unrestricted heteroskedasticity",
            "fr": "hétéroscédasticité non restreinte",
        },
        "diag.white.rows": {
            "en": "regressors must have one row per residual ({n}), got {got}.",
            "fr": "regressors doit avoir une ligne par résidu ({n}), reçu {got}.",
        },
        "diag.white.no_regressors": {
            "en": "White's test needs at least one regressor besides the constant.",
            "fr": "le test de White exige au moins un régresseur en plus de la constante.",
        },
        "diag.white.no_variation": {
            "en": "the regressors and their products are constant or collinear with the "
            "constant: nothing to test.",
            "fr": "les régresseurs et leurs produits sont constants ou colinéaires à la "
            "constante : rien à tester.",
        },
    },
)


@register_test(
    "white",
    TestSpec(
        family="heteroskedasticity",
        category="B",
        wave=1,
        reference="White, H. (1980), Econometrica 48(4), 817-838",
        stata="estat imtest",
        validation=("white_one", "white_three", "white_four", "white_homoskedastic"),
        inputs="arrays",
        balanced="na",
        r="lmtest::bptest(studentize=TRUE) with squares and products",
        summary="White's general test that the variance of the errors is constant.",
    ),
)
def white(residuals: object, regressors: object) -> TestResult:
    """Return White's test: ``N * R^2`` of ``e^2`` on the regressors, squares and products.

    The squared residuals are regressed on a constant, the regressors and all their squares and
    cross-products (centred, as Stata does for numerical stability); the statistic is ``N * R^2``,
    a chi-square with as many degrees of freedom as independent columns remain (constant or
    collinear products, such as the square of a dummy, are dropped).

    Parameters
    ----------
    residuals
        ``(n,)`` residuals of the fitted model.
    regressors
        ``(n,)`` or ``(n, p)`` regressors of the model, without the constant.

    Raises
    ------
    SpecificationError
        If the shapes do not match or there is no regressor.
    NumericalError
        If the auxiliary regressors carry no information.
    """
    # Stata: estat imtest, white -> regress e^2 on x and x_i * x_j (centred); chi2 = e(N) * e(r2);
    # df = e(df_m)
    e = as_array("residuals", residuals, ndim=1)
    raw = np.asarray(regressors, dtype=np.float64)
    x = as_array("regressors", raw[:, None] if raw.ndim == 1 else raw, ndim=2)
    if x.shape[0] != e.size:
        raise SpecificationError(
            _i18n.message("diag.white.rows", n=e.size, got=x.shape[0]), key="diag.white.rows"
        )
    if x.shape[1] == 0:
        raise SpecificationError(
            _i18n.message("diag.white.no_regressors"), key="diag.white.no_regressors"
        )
    centred = x - x.mean(axis=0)
    i, j = np.tril_indices(x.shape[1])
    auxiliary = np.column_stack([x, centred[:, i] * centred[:, j]])
    fit = ols_fit(e**2, auxiliary)
    if fit.df_model < 1:
        raise NumericalError(
            _i18n.message("diag.white.no_variation"), key="diag.white.no_variation"
        )
    statistic = e.size * fit.r2
    return TestResult(
        name="white",
        statistic=statistic,
        pvalue=float(stats.chi2.sf(statistic, fit.df_model)),
        statistic_name="chi2",
        df=(float(fit.df_model),),
        n_obs=e.size,
        options={"regressors": x.shape[1]},
    )
