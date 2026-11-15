"""Jarque-Bera test for normality, from bare residuals (category B, L.8.6).

Validated against Stata 17: Stata has no official command, so the case computes the statistic
in Stata from the skewness and kurtosis of ``summarize, detail`` and the published formula
(``n/6 * (S^2 + (K - 3)^2 / 4)``) and its chi-square tail, on residuals of several sizes; the
residuals are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_normality_jarque_bera.py``.
"""

from __future__ import annotations

from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, as_array, register_test
from pyxtpanel.diagnostics.normality.skewness_kurtosis import moments

_i18n.register(
    "diag",
    {
        "diag.jarque_bera.title": {
            "en": "Jarque-Bera test for normality",
            "fr": "Test de normalité de Jarque-Bera",
        },
        "diag.jarque_bera.null": {
            "en": "skewness is 0 and kurtosis is 3 (normality)",
            "fr": "l'asymétrie vaut 0 et l'aplatissement 3 (normalité)",
        },
        "diag.jarque_bera.alternative": {
            "en": "the data are not normally distributed",
            "fr": "les données ne suivent pas une loi normale",
        },
        "diag.jarque_bera.too_few": {
            "en": "n must be at least 2, got {n}.",
            "fr": "n doit valoir au moins 2, reçu {n}.",
        },
    },
)


@register_test(
    "jarque_bera",
    TestSpec(
        family="normality",
        category="B",
        wave=1,
        reference="Jarque, C. M. and Bera, A. K. (1987), International Statistical Review 55(2)",
        stata="summarize",
        validation=("jb_heavy_tails", "jb_normal", "jb_small"),
        inputs="arrays",
        balanced="na",
        r="tseries::jarque.bera.test",
        summary="Asymptotic joint test of zero skewness and kurtosis three.",
    ),
)
def jarque_bera(residuals: object) -> TestResult:
    """Return ``JB = n / 6 * (S^2 + (K - 3)^2 / 4)`` against a chi-square with 2 degrees of freedom.

    ``S`` and ``K`` are the skewness and the kurtosis with moments divided by ``n``.  The
    distribution is asymptotic: in small samples prefer :func:`skewness_kurtosis`, whose
    p-value is adjusted.

    Parameters
    ----------
    residuals
        ``(n,)`` values to test.

    Raises
    ------
    SpecificationError
        If there are fewer than 2 values.
    NumericalError
        If the values are constant.
    """
    # Stata: no official command; computed from summarize, detail (skewness, kurtosis, N)
    x = as_array("residuals", residuals, ndim=1)
    if x.size < 2:
        raise SpecificationError(
            _i18n.message("diag.jarque_bera.too_few", n=x.size), key="diag.jarque_bera.too_few"
        )
    skew, kurt = moments(x)
    statistic = x.size / 6.0 * (skew**2 + (kurt - 3.0) ** 2 / 4.0)
    return TestResult(
        name="jarque_bera",
        statistic=statistic,
        pvalue=float(stats.chi2.sf(statistic, 2)),
        statistic_name="chi2",
        df=(2.0,),
        n_obs=x.size,
        options={"skewness": skew, "kurtosis": kurt},
    )
