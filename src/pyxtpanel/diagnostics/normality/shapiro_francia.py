"""Shapiro-Francia W' test for normality, from bare residuals (category B, L.8.7).

Validated against Stata 17: ``sfrancia`` (W', V', z and the p-value, with the log transformation
of Royston 1993 and with ``boxcox``) on residuals of several sizes, with heavy tails, with
normal data and with many tied values; the residuals are exported from Stata and given to the
function.  The cases are in ``validation/cases/diagnostics_normality_shapiro_francia.py``.
"""

from __future__ import annotations

import math
from typing import Literal

import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, as_array, register_test
from pyxtpanel.diagnostics._shared.normal_scores import sorted_scores, squared_correlation

_i18n.register(
    "diag",
    {
        "diag.shapiro_francia.title": {
            "en": "Shapiro-Francia W' test for normality",
            "fr": "Test W' de normalité de Shapiro-Francia",
        },
        "diag.shapiro_francia.null": {
            "en": "the data are normally distributed",
            "fr": "les données suivent une loi normale",
        },
        "diag.shapiro_francia.alternative": {
            "en": "the data are not normally distributed",
            "fr": "les données ne suivent pas une loi normale",
        },
        "diag.shapiro_francia.too_few": {
            "en": "{n} observations: with the {transform} transformation the test needs at "
            "least {need}, the statistic is not computed.",
            "fr": "{n} observations : avec la transformation {transform} le test en exige au "
            "moins {need}, la statistique n'est pas calculée.",
        },
        "diag.shapiro_francia.range": {
            "en": "the normal approximation of W' is valid for {low} <= n <= {high} under the "
            "{transform} transformation; n = {n}.",
            "fr": "l'approximation normale de W' est valable pour {low} <= n <= {high} avec la "
            "transformation {transform} ; n = {n}.",
        },
        "diag.shapiro_francia.constant": {
            "en": "the data are constant: W' is not defined.",
            "fr": "les données sont constantes : W' n'est pas défini.",
        },
        "diag.shapiro_francia.transform": {
            "en": "transform must be 'log' or 'boxcox', got {value!r}.",
            "fr": "transform doit valoir 'log' ou 'boxcox', reçu {value!r}.",
        },
    },
)

#: transformation -> (fewest observations, largest valid n), as ``sfrancia``.
LIMITS = {"log": (10, 5000), "boxcox": (5, 1000)}


@register_test(
    "shapiro_francia",
    TestSpec(
        family="normality",
        category="B",
        wave=2,
        reference="Shapiro, S. S. and Francia, R. S. (1972), JASA 67(337), 215-216; "
        "Royston, P. (1983), Applied Statistics 32(2), 121-133; "
        "Royston, P. (1993), Statistician 42(1), 37-43",
        stata="sfrancia",
        validation=(
            "sf_heavy_tails",
            "sf_normal",
            "sf_boxcox",
            "sf_boxcox_small",
            "sf_small",
            "sf_ties",
            "sf_noties",
            "sf_single_precision",
        ),
        inputs="arrays",
        balanced="na",
        r="nortest::sf.test",
        aliases=("sfrancia",),
        summary="Squared correlation between the data and their normal scores.",
    ),
)
def shapiro_francia(
    residuals: object, *, ties: bool = True, transform: Literal["log", "boxcox"] = "log"
) -> TestResult:
    """Return Shapiro and Francia's ``W'`` with Royston's normal approximation.

    ``W'`` is the squared correlation between the data and the normal scores
    ``invnorm((rank - 0.375) / (n + 0.25))``; ``z`` is a normal transformation of it and the
    p-value its upper tail.  ``transform="log"`` (the default, as ``sfrancia`` since Stata 12) is
    Royston's (1993) approximation, valid for ``10 <= n <= 5000``; ``"boxcox"`` is Royston's
    (1983), valid for ``5 <= n <= 1000``.  ``details`` gives ``W'``, ``V'`` and ``z``.

    Parameters
    ----------
    residuals
        ``(n,)`` values to test.
    ties
        Average rank for tied values (the default of Stata); ``False`` is ``noties``.
    transform
        ``"log"`` or ``"boxcox"`` (Stata's ``boxcox`` option).

    Returns
    -------
    TestResult
        ``NaN`` with a remark when ``n`` is below the smallest size of the transformation; a
        remark when ``n`` exceeds the largest.

    Raises
    ------
    SpecificationError
        If ``transform`` is neither ``"log"`` nor ``"boxcox"``.
    NumericalError
        If the values are constant.
    """
    # Stata: sfrancia (sfrancia.ado): W' = corr(x, invnorm((rank - .375) / (n + .25)))^2, then the
    # log transformation of Royston (1993) or, with boxcox, the Box-Cox one of Royston (1983)
    if transform not in LIMITS:
        raise SpecificationError(
            _i18n.message("diag.shapiro_francia.transform", value=transform),
            key="diag.shapiro_francia.transform",
        )
    x = as_array("residuals", residuals, ndim=1)
    n = x.size
    low, high = LIMITS[transform]
    if n < low:
        return TestResult(
            name="shapiro_francia",
            statistic=math.nan,
            pvalue=math.nan,
            statistic_name="W'",
            n_obs=n,
            notes=(
                Note("diag.shapiro_francia.too_few", {"n": n, "transform": transform, "need": low}),
            ),
        )
    ordered, scores = sorted_scores(x, ties=ties)
    if ordered[0] == ordered[-1]:
        raise NumericalError(
            _i18n.message("diag.shapiro_francia.constant"), key="diag.shapiro_francia.constant"
        )
    r = squared_correlation(ordered, scores)
    v, z = _log(r, n) if transform == "log" else _boxcox(r, n)
    notes: tuple[Note, ...] = ()
    if n > high:
        notes = (
            Note(
                "diag.shapiro_francia.range",
                {"low": low, "high": high, "transform": transform, "n": n},
            ),
        )
    details = pd.DataFrame({"value": [r, v, z]}, index=pd.Index(["W'", "V'", "z"], name="quantity"))
    return TestResult(
        name="shapiro_francia",
        statistic=r,
        pvalue=float(stats.norm.sf(z)),
        statistic_name="W'",
        n_obs=n,
        options={"ties": ties, "transform": transform},
        details=details,
        notes=notes,
    )


def _log(r: float, n: int) -> tuple[float, float]:
    """Return ``V'`` and ``z`` under Royston's (1993) log transformation."""
    u = math.log(n)
    v = math.log(u)
    m = -1.2725 + 1.0521 * (v - u)
    s = 1.0308 - 0.26758 * (v + 2 / u)
    return (1 - r) * math.exp(-m), (math.log1p(-r) - m) / s


def _boxcox(r: float, n: int) -> tuple[float, float]:
    """Return ``V'`` and ``z`` under Royston's (1983) Box-Cox transformation."""
    x = math.log(n) - 5
    lam = -0.0480157 + x * (0.01971964 - 0.0119065 * x * x)
    m = -math.exp(
        1.6930674 + x * (0.1441647 + x * (-0.01849276 + x * (0.031074485 + x * 0.0055717663)))
    )
    s = math.exp(
        -0.510725 + x * (-0.1160364 + x * (-0.006702098 + x * (0.054465944 + x * 0.0087397329)))
    )
    y = ((1 - r) ** lam - 1) / lam
    return (1 - r) / (lam * m + 1) ** (1 / lam), (y - m) / s
