"""Sargan test of over-identifying restrictions, from bare arrays (category B, L.6.2).

Validated against Stata 17: ``ivregress 2sls`` followed by ``estat overid`` (Sargan) and
``xtabond2`` (``e(sargan)``); the cases are in
``validation/cases/diagnostics_postestimation_sargan.py``.
"""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.diagnostics._core import TestResult, TestSpec, register_test
from pyxtpanel.diagnostics._shared.overid import overidentification

_i18n.register(
    "diag",
    {
        "diag.sargan.title": {
            "en": "Sargan test of overidentifying restrictions",
            "fr": "Test de Sargan de restrictions de sur-identification",
        },
        "diag.sargan.null": {
            "en": "the overidentifying restrictions are valid",
            "fr": "les restrictions de sur-identification sont valides",
        },
        "diag.sargan.alternative": {
            "en": "some overidentifying restrictions are not valid",
            "fr": "certaines restrictions de sur-identification ne sont pas valides",
        },
    },
)


@register_test(
    "sargan",
    TestSpec(
        family="postestimation",
        category="B",
        wave=1,
        reference="Sargan, J. D. (1958), Econometrica 26(3), 393-415",
        stata="ivregress 2sls; xtabond2",
        validation=(
            "sargan_ivregress",
            "sargan_ivregress_invalid",
            "sargan_just_identified",
            "sargan_ab1_balanced",
            "sargan_ab2_balanced",
            "sargan_ab2_gaps",
            "sargan_ab2_entry_exit",
        ),
        inputs="arrays",
        balanced="na",
        r="plm::sargan",
        summary="Sargan test of the overidentifying restrictions of a one-step fit.",
    ),
)
def sargan(
    moments: object,
    weight: object,
    n_parameters: int,
    *,
    scale: float,
    n_instruments: int | None = None,
    n_groups: int | None = None,
) -> TestResult:
    """Return Sargan's statistic, ``m' W m / scale`` against a chi-square with ``q - k`` df.

    ``m = sum_i Z_i' e_i`` are the moments at the one-step estimate, ``W`` the weight matrix of
    that step (``(Z'Z)^-1`` for two-stage least squares) and ``scale`` the residual variance the
    estimator uses (``e'e / N`` for ``ivregress``; ``xtabond2`` has its own convention for the
    differenced residuals).  The estimator builds all three.

    Parameters
    ----------
    moments
        ``(q,)`` vector ``sum_i Z_i' e_i``, or the ``(n_units, q)`` moments of each unit.
    weight
        ``(q, q)`` weight matrix of the one-step criterion (symmetrised).
    n_parameters
        Number of estimated parameters.
    scale
        Positive residual variance dividing the quadratic form.
    n_instruments
        Independent instruments, ``q`` by default (pass the rank when some are collinear).
    n_groups
        Number of units: a remark is added when instruments outnumber them.

    Returns
    -------
    TestResult
        ``statistic`` and ``pvalue`` are ``NaN`` when the model is not over-identified.
    """
    # Stata: ivregress 2sls ... + estat overid -> r(sargan), r(df), r(p_sargan)
    # Stata: xtabond2 -> e(sargan), e(sarganp)
    return overidentification(
        "sargan",
        moments,
        weight,
        n_parameters,
        scale=scale,
        n_instruments=n_instruments,
        n_groups=n_groups,
    )
