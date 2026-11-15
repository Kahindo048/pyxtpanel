"""Difference-in-Hansen (C) test of the exogeneity of a subset of instruments (category B, L.6.3).

Validated against Stata 17: the "Difference-in-Hansen tests of exogeneity of instrument subsets"
of ``xtabond2`` (``e(diffsargan)``), for one and for several instrument groups, on the balanced
panel and on a panel with gaps; the cases are in
``validation/cases/diagnostics_postestimation_diff_hansen.py``.
"""

from __future__ import annotations

import math

from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, register_test

_i18n.register(
    "diag",
    {
        "diag.diff_hansen.title": {
            "en": "Difference-in-Hansen test of exogeneity of instrument subsets",
            "fr": "Test de Hansen en différence d'exogénéité d'un sous-ensemble d'instruments",
        },
        "diag.diff_hansen.null": {
            "en": "the instruments of the subset are exogenous",
            "fr": "les instruments du sous-ensemble sont exogènes",
        },
        "diag.diff_hansen.alternative": {
            "en": "some instruments of the subset are not exogenous",
            "fr": "certains instruments du sous-ensemble ne sont pas exogènes",
        },
        "diag.diff_hansen.negative": {
            "en": "the difference is negative ({value:.4g}): the restricted statistic exceeds the "
            "full one, so the subset test is not valid for these fits; no p-value.",
            "fr": "la différence est négative ({value:.4g}) : la statistique restreinte dépasse "
            "la complète, le test du sous-ensemble n'est pas valide pour ces estimations ; "
            "pas de p-value.",
        },
        "diag.diff_hansen.bad_argument": {
            "en": "{name} must be {expected}, got {value}.",
            "fr": "{name} doit être {expected}, reçu {value}.",
        },
    },
)


@register_test(
    "diff_hansen",
    TestSpec(
        family="postestimation",
        category="B",
        wave=2,
        reference="Eichenbaum, M. S., Hansen, L. P. and Singleton, K. J. (1988), QJE 103(1), 51-78",
        stata="xtabond2",
        validation=("diffhansen_balanced", "diffhansen_gaps_groups"),
        inputs="arrays",
        balanced="na",
        r=None,
        summary="Difference of two Hansen statistics: exogeneity of a subset of instruments.",
    ),
)
def diff_hansen(
    j_full: float,
    j_excluding: float,
    n_dropped: int,
    *,
    df_full: int | None = None,
) -> TestResult:
    """Return the C statistic ``J_full - J_excluding`` against a chi-square with ``n_dropped`` df.

    ``J_full`` is the Hansen statistic with all the instruments and ``J_excluding`` the one
    computed without the subset under test, with the weight matrix of the full model restricted
    to the remaining instruments (as ``xtabond2`` does, so that the difference is not negative in
    theory).  The estimator computes both; this function takes the difference and its
    distribution.  A negative difference, possible in finite samples, is returned as such with no
    p-value (the package's convention; Stata's outputs checked here are all non-negative).

    Parameters
    ----------
    j_full, j_excluding
        The two Hansen statistics, non-negative.
    n_dropped
        Number of instruments in the subset (the degrees of freedom of the difference).
    df_full
        Degrees of freedom of ``j_full``, to check that the restricted model stays identified.

    Raises
    ------
    SpecificationError
        If a statistic is negative or not finite, ``n_dropped`` is not a positive integer, or
        the subset is larger than the over-identification of the full model.
    """
    # Stata: xtabond2 -> e(diffsargan): "Hansen test excluding group", "Difference (null H =
    # exogenous)" = J_full - J_excluding on as many degrees of freedom as instruments dropped
    for name, value in (("j_full", j_full), ("j_excluding", j_excluding)):
        if not (math.isfinite(value) and value >= 0):
            raise SpecificationError(
                _i18n.message(
                    "diag.diff_hansen.bad_argument",
                    name=name,
                    expected="a finite non-negative statistic",
                    value=value,
                ),
                key="diag.diff_hansen.bad_argument",
            )
    if float(n_dropped) != int(n_dropped) or n_dropped < 1:
        raise SpecificationError(
            _i18n.message(
                "diag.diff_hansen.bad_argument",
                name="n_dropped",
                expected="a positive integer",
                value=n_dropped,
            ),
            key="diag.diff_hansen.bad_argument",
        )
    if df_full is not None and n_dropped > df_full:
        raise SpecificationError(
            _i18n.message(
                "diag.diff_hansen.bad_argument",
                name="n_dropped",
                expected=f"at most the {df_full} degrees of freedom of the full model",
                value=n_dropped,
            ),
            key="diag.diff_hansen.bad_argument",
        )
    value = j_full - j_excluding
    notes: list[Note] = []
    if value < 0:
        notes.append(Note("diag.diff_hansen.negative", {"value": value}))
        pvalue = math.nan
    else:
        pvalue = float(stats.chi2.sf(value, n_dropped))
    return TestResult(
        name="diff_hansen",
        statistic=value,
        pvalue=pvalue,
        statistic_name="chi2",
        df=(float(n_dropped),),
        options={"j_full": j_full, "j_excluding": j_excluding},
        notes=tuple(notes),
    )
