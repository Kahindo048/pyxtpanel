"""Skewness and kurtosis test for normality, from bare residuals (category B, L.8.6).

Validated against Stata 17: ``sktest`` (the adjusted and the unadjusted joint chi-square, and
the two separate p-values) on residuals of several sizes, including 9 observations; the
residuals are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_normality_skewness_kurtosis.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, as_array, register_test

_i18n.register(
    "diag",
    {
        "diag.skewness_kurtosis.title": {
            "en": "Skewness and kurtosis test for normality",
            "fr": "Test d'asymétrie et d'aplatissement de normalité",
        },
        "diag.skewness_kurtosis.null": {
            "en": "the data are normally distributed",
            "fr": "les données suivent une loi normale",
        },
        "diag.skewness_kurtosis.alternative": {
            "en": "the data are not normally distributed",
            "fr": "les données ne suivent pas une loi normale",
        },
        "diag.skewness_kurtosis.too_few": {
            "en": "{n} observations: the test needs at least 8, the statistic is not computed.",
            "fr": "{n} observations : le test en exige au moins 8, la statistique n'est pas "
            "calculée.",
        },
        "diag.skewness_kurtosis.undefined": {
            "en": "the kurtosis is so far below the normal value that its transformation is not "
            "defined (Stata reports a missing value too).",
            "fr": "l'aplatissement est si inférieur à la valeur normale que sa transformation "
            "n'est pas définie (Stata renvoie aussi une valeur manquante).",
        },
        "diag.skewness_kurtosis.constant": {
            "en": "the residuals are constant: skewness and kurtosis are not defined.",
            "fr": "les résidus sont constants : l'asymétrie et l'aplatissement ne sont pas "
            "définis.",
        },
    },
)

#: Fewest observations for which the test is computed (as ``sktest``).
MIN_OBS = 8


def moments(x: np.ndarray) -> tuple[float, float]:
    """Return the skewness ``m3 / m2^1.5`` and the kurtosis ``m4 / m2^2`` (divisor ``n``).

    Raises
    ------
    NumericalError
        If ``x`` is constant.
    """
    centred = x - x.mean()
    m2 = float(np.mean(centred**2))
    if m2 == 0.0:
        raise NumericalError(
            _i18n.message("diag.skewness_kurtosis.constant"),
            key="diag.skewness_kurtosis.constant",
        )
    return float(np.mean(centred**3)) / m2**1.5, float(np.mean(centred**4)) / m2**2


@register_test(
    "skewness_kurtosis",
    TestSpec(
        family="normality",
        category="B",
        wave=1,
        reference="D'Agostino, R. B., Belanger, A. and D'Agostino Jr., R. B. (1990), "
        "The American Statistician 44(4), 316-321; Royston, P. (1991), Stata Technical Bulletin 3",
        stata="sktest",
        validation=("sk_heavy_tails", "sk_normal", "sk_unadjusted", "sk_small"),
        inputs="arrays",
        balanced="na",
        r="moments::agostino.test / anscombe.test",
        summary="Joint test of normality from the skewness and the kurtosis of the residuals.",
    ),
)
def skewness_kurtosis(residuals: object, *, adjust: bool = True) -> TestResult:
    """Return the D'Agostino-Belanger-D'Agostino joint test of skewness and kurtosis.

    The skewness and the kurtosis are turned into two standard normal variables ``Z1`` and
    ``Z2`` (D'Agostino's transformation of the skewness, Anscombe-Glynn's of the kurtosis) and
    the statistic is ``Z1^2 + Z2^2``, a chi-square with 2 degrees of freedom.  With ``adjust``
    (the default, as ``sktest``) Royston's empirical adjustment is applied so that the p-value
    matches the tables.  ``details`` gives the skewness, the kurtosis and their separate tests.

    Parameters
    ----------
    residuals
        ``(n,)`` values to test, ``n >= 8``.
    adjust
        Royston's adjustment (``sktest`` without ``noadjust``).

    Returns
    -------
    TestResult
        ``NaN`` with a remark when ``n < 8``.
    """
    # Stata: sktest -> Mata SKTest (sktest.ado): Z1 (skewness), Z2 (kurtosis), K2 = Z1^2 + Z2^2,
    # Royston's adjustment unless noadjust; moments from summarize, detail (divisor n)
    x = as_array("residuals", residuals, ndim=1)
    n = x.size
    if n < MIN_OBS:
        return TestResult(
            name="skewness_kurtosis",
            statistic=math.nan,
            pvalue=math.nan,
            statistic_name="chi2",
            df=(2.0,),
            n_obs=n,
            notes=(Note("diag.skewness_kurtosis.too_few", {"n": n}),),
        )
    skew, kurt = moments(x)
    z1 = _z_skewness(skew, n)
    z2 = _z_kurtosis(kurt, n)
    k2 = z1 * z1 + z2 * z2
    if not math.isfinite(k2):
        return TestResult(
            name="skewness_kurtosis",
            statistic=math.nan,
            pvalue=math.nan,
            statistic_name="chi2",
            df=(2.0,),
            n_obs=n,
            notes=(Note("diag.skewness_kurtosis.undefined"),),
        )
    if adjust:
        k2, pvalue = _royston(k2, n)
    else:
        pvalue = float(stats.chi2.sf(k2, 2))
    details = pd.DataFrame(
        {
            "value": [skew, kurt],
            "z": [z1, z2],
            "pvalue": [2 * float(stats.norm.sf(abs(z1))), 2 * float(stats.norm.sf(abs(z2)))],
        },
        index=pd.Index(["skewness", "kurtosis"], name="moment"),
    )
    return TestResult(
        name="skewness_kurtosis",
        statistic=k2,
        pvalue=pvalue,
        statistic_name="chi2",
        df=(2.0,),
        n_obs=n,
        options={"adjust": adjust},
        details=details,
    )


def _z_skewness(skewness: float, n: int) -> float:
    y = skewness * math.sqrt(((n + 1) * (n + 3)) / (6 * (n - 2)))
    beta2 = (3 * (n * n + 27 * n - 70) * (n + 1) * (n + 3)) / (
        (n - 2) * (n + 5) * (n + 7) * (n + 9)
    )
    w2 = -1 + math.sqrt(2 * (beta2 - 1))
    delta = 1 / math.sqrt(math.log(math.sqrt(w2)))
    alpha = math.sqrt(2 / (w2 - 1))
    return delta * math.log(y / alpha + math.sqrt((y / alpha) ** 2 + 1))


def _z_kurtosis(kurtosis: float, n: int) -> float:
    eb2 = (3 * (n - 1)) / (n + 1)
    vb2 = (24 * n * (n - 2) * (n - 3)) / ((n + 1) ** 2 * (n + 3) * (n + 5))
    x = (kurtosis - eb2) / math.sqrt(vb2)
    rbeta1 = ((6 * (n * n - 5 * n + 2)) / ((n + 7) * (n + 9))) * math.sqrt(
        (6 * (n + 3) * (n + 5)) / (n * (n - 2) * (n - 3))
    )
    a = 6 + (8 / rbeta1) * (2 / rbeta1 + math.sqrt(1 + 4 / (rbeta1 * rbeta1)))
    cube = (1 - 2 / a) / (1 + x * math.sqrt(2 / (a - 4)))
    if cube < 0:  # Mata's negative ^ (1/3) is missing: sktest reports nothing either
        return math.nan
    return ((1 - 2 / (9 * a)) - cube ** (1 / 3)) / math.sqrt(2 / (9 * a))


def _royston(k2: float, n: int) -> tuple[float, float]:
    """Royston's empirical adjustment of the chi-square(2) statistic and its p-value."""
    zc2 = -float(stats.norm.ppf(math.exp(-0.5 * k2)))
    logn = math.log(n)
    cut = 0.55 * n**0.2 - 0.21
    a1 = (-5 + 3.46 * logn) * math.exp(-1.37 * logn)
    b1 = 1 + (0.854 - 0.148 * logn) * math.exp(-0.55 * logn)
    b2mb1 = 2.13 / (1 - 2.37 * logn)
    a2 = a1 - b2mb1 * cut
    b2 = b2mb1 + b1
    if zc2 < -1:
        z = zc2
    elif zc2 < cut:
        z = a1 + b1 * zc2
    else:
        z = a2 + b2 * zc2
    return -2.0 * float(stats.norm.logcdf(-z)), float(stats.norm.sf(z))
