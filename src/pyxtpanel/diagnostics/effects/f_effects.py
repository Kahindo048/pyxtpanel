"""F test that the unit effects are all zero, from two sums of squares (category B, L.8.5).

Validated against Stata 17: the ``F test that all u_i=0`` of ``xtreg, fe`` (``e(F_f)``), with
the pooled residual sum of squares taken from ``regress``, on a balanced panel, an unbalanced
panel and a panel with gaps.  The cases are in
``validation/cases/diagnostics_effects_f_effects.py``.
"""

from __future__ import annotations

import math

from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import TestResult, TestSpec, register_test

_i18n.register(
    "diag",
    {
        "diag.f_effects.title": {
            "en": "F test that all unit effects are zero",
            "fr": "Test F de nullité de tous les effets individuels",
        },
        "diag.f_effects.null": {
            "en": "u_i = 0 for all units i",
            "fr": "u_i = 0 pour toute unité i",
        },
        "diag.f_effects.alternative": {
            "en": "at least one unit effect is not zero",
            "fr": "au moins un effet individuel n'est pas nul",
        },
        "diag.f_effects.bad_argument": {
            "en": "{name} must be {expected}, got {value}.",
            "fr": "{name} doit être {expected}, reçu {value}.",
        },
    },
)


@register_test(
    "f_effects",
    TestSpec(
        family="effects",
        category="B",
        wave=1,
        reference="Baltagi, B. H. (2021), Econometric Analysis of Panel Data, 6th ed., ch. 4",
        stata="xtreg",
        validation=("fe_balanced", "fe_unbalanced", "fe_gaps"),
        inputs="arrays",
        balanced="na",
        r="plm::pFtest",
        summary="F test of the fixed effects from the pooled and within residual sums of squares.",
    ),
)
def f_effects(
    rss_pooled: float,
    rss_effects: float,
    n_obs: int,
    n_groups: int,
    n_regressors: int,
) -> TestResult:
    """Return ``F = ((RSS_pooled - RSS_fe) / (G - 1)) / (RSS_fe / (N - G - k))``.

    ``RSS_pooled`` is the residual sum of squares of the pooled regression on a constant and the
    ``k`` regressors, ``RSS_fe`` the one of the within (fixed-effects) regression.  The statistic
    follows an ``F(G - 1, N - G - k)`` under the null that no unit has its own intercept.

    Parameters
    ----------
    rss_pooled, rss_effects
        The two residual sums of squares, ``rss_pooled >= rss_effects > 0``.
    n_obs, n_groups
        Number of observations ``N`` and of units ``G``.
    n_regressors
        Number ``k`` of slope coefficients (without the constant).

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If the sums of squares are inconsistent or the degrees of freedom are not positive.
    """
    # Stata: xtreg, fe -> e(F_f) with e(df_a) = G - 1 and e(df_r) = N - G - k
    for name, value in (("rss_pooled", rss_pooled), ("rss_effects", rss_effects)):
        if not (math.isfinite(value) and value > 0):
            raise _error(name, "a positive finite number", value)
    if rss_pooled < rss_effects:
        raise _error("rss_pooled", f"at least rss_effects ({rss_effects})", rss_pooled)
    df1, df2 = n_groups - 1, n_obs - n_groups - n_regressors
    if df1 < 1:
        raise _error("n_groups", "at least 2", n_groups)
    if df2 < 1:
        raise _error(
            "n_obs", f"larger than n_groups + n_regressors ({n_groups + n_regressors})", n_obs
        )
    statistic = ((rss_pooled - rss_effects) / df1) / (rss_effects / df2)
    return TestResult(
        name="f_effects",
        statistic=statistic,
        pvalue=float(stats.f.sf(statistic, df1, df2)),
        statistic_name="F",
        df=(float(df1), float(df2)),
        n_obs=n_obs,
        n_groups=n_groups,
        options={"n_regressors": n_regressors},
    )


def _error(name: str, expected: str, value: object) -> SpecificationError:
    return SpecificationError(
        _i18n.message("diag.f_effects.bad_argument", name=name, expected=expected, value=value),
        key="diag.f_effects.bad_argument",
    )
