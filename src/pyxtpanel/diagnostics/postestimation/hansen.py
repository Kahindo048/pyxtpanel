"""Hansen J test of over-identifying restrictions, from bare arrays (category B, L.6.2).

Validated against Stata 17: ``ivregress gmm ..., wmatrix(cluster id)`` followed by
``estat overid`` (Hansen's J), ``xtabond2 ..., twostep`` (``e(hansen)``) and ``xtdpthresh``;
the cases are in ``validation/cases/diagnostics_postestimation_hansen.py``.
"""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.diagnostics._core import TestResult, TestSpec, register_test
from pyxtpanel.diagnostics._shared.overid import overidentification

_i18n.register(
    "diag",
    {
        "diag.hansen.title": {
            "en": "Hansen J test of overidentifying restrictions",
            "fr": "Test J de Hansen de restrictions de sur-identification",
        },
        "diag.hansen.null": {
            "en": "the overidentifying restrictions are valid",
            "fr": "les restrictions de sur-identification sont valides",
        },
        "diag.hansen.alternative": {
            "en": "some overidentifying restrictions are not valid",
            "fr": "certaines restrictions de sur-identification ne sont pas valides",
        },
    },
)


@register_test(
    "hansen",
    TestSpec(
        family="postestimation",
        category="B",
        wave=1,
        reference="Hansen, L. P. (1982), Econometrica 50(4), 1029-1054",
        stata="ivregress gmm; xtabond2",
        validation=(
            "hansen_ivregress",
            "hansen_ivregress_invalid",
            "hansen_just_identified",
            "hansen_ab2_balanced",
            "hansen_ab2_gaps",
            "hansen_ab2_entry_exit",
        ),
        inputs="arrays",
        balanced="na",
        r=None,
        summary="J test of the overidentifying restrictions of an efficient GMM fit.",
    ),
)
def hansen(
    moments: object,
    weight: object,
    n_parameters: int,
    *,
    n_instruments: int | None = None,
    n_groups: int | None = None,
) -> TestResult:
    """Return Hansen's J test, ``m' W m`` against a chi-square with ``q - k`` degrees of freedom.

    ``m = sum_i Z_i' e_i`` are the moments at the final estimate and ``W`` the weight matrix of
    the criterion, the inverse of the covariance of the moments (``moment_covariance`` gives
    ``sum_i g_i g_i'`` from the moments of each unit).  The estimator builds both: which
    residuals, which instruments and which weight matrix are conventions of the estimator.

    Parameters
    ----------
    moments
        ``(q,)`` vector ``sum_i Z_i' e_i``, or the ``(n_units, q)`` moments of each unit.
    weight
        ``(q, q)`` weight matrix of the criterion (symmetrised).
    n_parameters
        Number of estimated parameters (count the threshold of a threshold model).
    n_instruments
        Independent instruments, ``q`` by default (pass the rank when some are collinear).
    n_groups
        Number of units: a remark is added when instruments outnumber them.

    Returns
    -------
    TestResult
        ``statistic`` and ``pvalue`` are ``NaN`` when the model is not over-identified.
    """
    # Stata: ivregress gmm ..., wmatrix(cluster id) + estat overid -> r(HansenJ), r(J_df)
    # Stata: xtabond2 ..., twostep -> e(hansen), e(hansen_df), e(hansenp)
    return overidentification(
        "hansen",
        moments,
        weight,
        n_parameters,
        n_instruments=n_instruments,
        n_groups=n_groups,
    )
