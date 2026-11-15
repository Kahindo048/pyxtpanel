"""Shapiro-Wilk W test for normality, from bare residuals (category B, L.8.7).

Validated against Stata 17: ``swilk`` (W, V, z and the p-value) on residuals of several sizes,
from 3 to 600 observations, with heavy tails, with normal data and with many tied values; the
residuals are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_normality_shapiro_wilk.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, as_array, register_test
from pyxtpanel.diagnostics._shared.normal_scores import sorted_scores, squared_correlation

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray

_i18n.register(
    "diag",
    {
        "diag.shapiro_wilk.title": {
            "en": "Shapiro-Wilk W test for normality",
            "fr": "Test W de normalité de Shapiro-Wilk",
        },
        "diag.shapiro_wilk.null": {
            "en": "the data are normally distributed",
            "fr": "les données suivent une loi normale",
        },
        "diag.shapiro_wilk.alternative": {
            "en": "the data are not normally distributed",
            "fr": "les données ne suivent pas une loi normale",
        },
        "diag.shapiro_wilk.too_few": {
            "en": "{n} observations: the test needs at least 3, the statistic is not computed.",
            "fr": "{n} observations : le test en exige au moins 3, la statistique n'est pas "
            "calculée.",
        },
        "diag.shapiro_wilk.range": {
            "en": "the normal approximation of W is valid for 4 <= n <= 2000; n = {n}.",
            "fr": "l'approximation normale de W est valable pour 4 <= n <= 2000 ; n = {n}.",
        },
        "diag.shapiro_wilk.constant": {
            "en": "the data are constant: W is not defined.",
            "fr": "les données sont constantes : W n'est pas défini.",
        },
        "diag.shapiro_wilk.boundary": {
            "en": "n = 3 and W at its lower bound: the approximate exact p-value is slightly "
            "negative and is reported as 0 (Stata prints it as it is and leaves z missing).",
            "fr": "n = 3 et W à sa borne inférieure : la p-value exacte approchée est "
            "légèrement négative et est ramenée à 0 (Stata l'affiche telle quelle et laisse z "
            "manquant).",
        },
    },
)

#: Stata's bounds of validity of the normal approximation of W.
MIN_VALID, MAX_VALID = 4, 2000


@register_test(
    "shapiro_wilk",
    TestSpec(
        family="normality",
        category="B",
        wave=2,
        reference="Shapiro, S. S. and Wilk, M. B. (1965), Biometrika 52(3-4), 591-611; "
        "Royston, P. (1982), Applied Statistics 31(2), 115-124 and 31(2), 176-180",
        stata="swilk",
        validation=(
            "sw_heavy_tails",
            "sw_normal",
            "sw_small",
            "sw_n3",
            "sw_n4",
            "sw_n5",
            "sw_n7",
            "sw_n11",
            "sw_n12",
            "sw_ties",
            "sw_noties",
            "sw_large",
            "sw_single_precision",
        ),
        inputs="arrays",
        balanced="na",
        r="stats::shapiro.test (Royston 1995: close, not identical)",
        aliases=("swilk",),
        summary="Squared correlation between the ordered data and the normal order statistics.",
    ),
)
def shapiro_wilk(residuals: object, *, ties: bool = True) -> TestResult:
    """Return Shapiro and Wilk's ``W`` with Royston's normal approximation of its distribution.

    ``W`` is the squared correlation between the ordered data and coefficients built from the
    expected normal order statistics (Royston's approximations of the coefficients, as ``swilk``
    computes them).  ``z`` is a normal transformation of ``log(1 - W)`` and the p-value its upper
    tail; ``V = (1 - W) / median(1 - W)`` is an index of departure from normality (1 when the
    data are normal).  For ``n = 3`` the p-value is exact.  ``details`` gives ``W``, ``V`` and
    ``z``.

    Parameters
    ----------
    residuals
        ``(n,)`` values to test, ``n >= 3``.
    ties
        Average rank for tied values (the default of Stata); ``False`` is ``swilk``'s
        ``noties`` option.

    Returns
    -------
    TestResult
        ``NaN`` with a remark when ``n < 3``; a remark when ``n`` is outside ``4..2000``.

    Raises
    ------
    NumericalError
        If the values are constant.
    """
    # Stata: swilk (swilk.ado): scores invnorm((rank - .375) / (n + .25)), Royston's polynomial
    # edge coefficients, W = corr^2, then (log(1 - W) - M) / S with Royston's M and S by n
    x = as_array("residuals", residuals, ndim=1)
    n = x.size
    if n < 3:
        return TestResult(
            name="shapiro_wilk",
            statistic=math.nan,
            pvalue=math.nan,
            statistic_name="W",
            n_obs=n,
            notes=(Note("diag.shapiro_wilk.too_few", {"n": n}),),
        )
    ordered, scores = sorted_scores(x, ties=ties)
    if ordered[0] == ordered[-1]:
        raise NumericalError(
            _i18n.message("diag.shapiro_wilk.constant"), key="diag.shapiro_wilk.constant"
        )
    w = squared_correlation(ordered, _coefficients(scores))
    notes: list[Note] = []
    if n < MIN_VALID or n > MAX_VALID:
        notes.append(Note("diag.shapiro_wilk.range", {"n": n}))
    if n == 3:
        v, z, pvalue = _exact_three(w)
        if pvalue < 0:
            notes.append(Note("diag.shapiro_wilk.boundary"))
            pvalue = 0.0
    else:
        v, z = _normalised(w, n)
        pvalue = float(stats.norm.sf(z))
    details = pd.DataFrame({"value": [w, v, z]}, index=pd.Index(["W", "V", "z"], name="quantity"))
    return TestResult(
        name="shapiro_wilk",
        statistic=w,
        pvalue=pvalue,
        statistic_name="W",
        n_obs=n,
        options={"ties": ties},
        details=details,
        notes=tuple(notes),
    )


def _coefficients(scores: FloatArray) -> FloatArray:
    """Return Royston's approximate coefficients ``a`` from the sorted normal scores."""
    n = scores.size
    if n == 3:
        return np.array([-math.sqrt(0.5), 0.0, math.sqrt(0.5)])
    summ2 = float(np.var(scores, ddof=1)) * (n - 1)
    x = 1 / math.sqrt(n)
    a = scores.copy()
    a1 = scores[-1] / math.sqrt(summ2) + x * (
        0.221157 + x * (-0.147981 + x * (-2.071190 + x * (4.434685 - x * 2.706056)))
    )
    if n > 5:
        a2 = scores[-2] / math.sqrt(summ2) + x * (
            0.042981 + x * (-0.293762 + x * (-1.752461 + x * (5.682633 - x * 3.582633)))
        )
        fac = math.sqrt(
            (summ2 - 2 * scores[-1] ** 2 - 2 * scores[-2] ** 2) / (1 - 2 * a1**2 - 2 * a2**2)
        )
        first = 2  # zero-based index of the first coefficient left to the plain scores
        a[-1], a[-2], a[0], a[1] = a1, a2, -a1, -a2
    else:
        fac = math.sqrt((summ2 - 2 * scores[-1] ** 2) / (1 - 2 * a1**2))
        first = 1
        a[-1], a[0] = a1, -a1
    a[first : n - first] = a[first : n - first] / fac
    return np.asarray(a, dtype=np.float64)


def _exact_three(w: float) -> tuple[float, float, float]:
    """Return ``V``, ``z`` and the p-value for ``n = 3`` (exact, Royston's approximation)."""
    sw = math.sqrt(w)
    ang = 1.5707288 + sw * (-0.2121144 + sw * (0.074261 - 0.0187293 * sw))
    ang = math.pi / 2 - ang * math.sqrt(1 - sw)
    stqr = 1.047198  # arcsin(sqrt(0.75))
    pvalue = (6 / math.pi) * (ang - stqr)
    z = -float(stats.norm.ppf(pvalue)) if 0 < pvalue < 1 else math.nan
    v = (1 - w) / (1 - math.sin(math.pi / 12 + stqr) ** 2)
    return v, z, pvalue


def _normalised(w: float, n: int) -> tuple[float, float]:
    """Return ``V`` and ``z`` for ``4 <= n`` from Royston's moments of ``log(1 - W)``."""
    y = math.log1p(-w)
    x = math.log(n)
    if n <= 11:
        gamma = n * 0.459 - 2.273
        if y >= gamma:
            return math.nan, 9.9999  # Stata's sentinel (M = 0, S = 1): z is a large constant
        y = -math.log(gamma - y)
        m = 0.5440 + n * (-0.39978 + n * (0.025054 - n * 0.0006714))
        s = math.exp(1.3822 + n * (-0.77857 + n * (0.062767 - n * 0.0020322)))
        v = (1 - w) / math.exp(gamma - math.exp(-m))
    else:
        m = -1.5861 + x * (-0.31082 + x * (-0.083751 + x * 0.0038915))
        s = math.exp(-0.4803 + x * (-0.082676 + x * 0.0030302))
        v = (1 - w) / math.exp(m)
    return v, (y - m) / s
