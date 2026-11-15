"""Breusch-Pagan Lagrange multiplier test for random effects, from bare arrays (L.8.4).

Validated against Stata 17: ``xttest0`` after ``xtreg, re`` on panels with and without a unit
effect (``sigma_u`` zero and positive), balanced and unbalanced; the pooled residuals and the
unit boundaries are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_effects_lm_random_effects.py``.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, as_array, register_test

_i18n.register(
    "diag",
    {
        "diag.lm_random_effects.title": {
            "en": "Breusch and Pagan Lagrangian multiplier test for random effects",
            "fr": "Test du multiplicateur de Lagrange de Breusch et Pagan d'effets aléatoires",
        },
        "diag.lm_random_effects.null": {"en": "Var(u) = 0", "fr": "Var(u) = 0"},
        "diag.lm_random_effects.alternative": {"en": "Var(u) > 0", "fr": "Var(u) > 0"},
        "diag.lm_random_effects.zero_effect": {
            "en": "the random-effects estimate of Var(u) is zero: the statistic is 0 and the "
            "p-value 1, as in Stata.",
            "fr": "l'estimation de Var(u) par les effets aléatoires est nulle : la statistique "
            "vaut 0 et la p-value 1, comme dans Stata.",
        },
        "diag.lm_random_effects.offsets": {
            "en": "offsets must start at 0, increase and end at the number of residuals ({n}).",
            "fr": "offsets doit commencer à 0, croître et finir au nombre de résidus ({n}).",
        },
        "diag.lm_random_effects.single": {
            "en": "every unit has a single observation: the test needs repeated observations.",
            "fr": "chaque unité n'a qu'une observation : le test exige des observations répétées.",
        },
        "diag.lm_random_effects.zero_residuals": {
            "en": "all residuals are zero: the test is not defined.",
            "fr": "tous les résidus sont nuls : le test n'est pas défini.",
        },
    },
)


@register_test(
    "lm_random_effects",
    TestSpec(
        family="effects",
        category="B",
        wave=1,
        reference="Breusch, T. S. and Pagan, A. R. (1980), Rev. Econ. Studies 47(1), 239-253",
        stata="xttest0",
        validation=("lm_re_effect", "lm_re_unbalanced", "lm_re_no_effect", "lm_re_gaps"),
        inputs="arrays",
        balanced="na",
        r="plm::plmtest(type='bp')",
        summary="Test that the variance of the unit effects is zero, from pooled residuals.",
    ),
)
def lm_random_effects(
    residuals: object, offsets: object, *, variance_u: float | None = None
) -> TestResult:
    """Return the Breusch-Pagan LM statistic ``N^2/2 * A^2 / (sum T_i^2 - N)`` for ``Var(u) = 0``.

    ``A = 1 - sum_i (sum_t e_it)^2 / sum e_it^2`` is built from the residuals of the pooled
    least-squares regression.  The p-value is one-sided (half of the chi-square(1) tail, the
    chibar-squared distribution).

    Parameters
    ----------
    residuals
        ``(n,)`` residuals of the pooled regression of the dependent variable on the regressors
        and a constant, sorted by unit.
    offsets
        ``(G + 1,)`` start of each unit in ``residuals``.
    variance_u
        The random-effects estimate of ``Var(u)``, when the estimator has it.  ``xttest0``
        reports a statistic of 0 and a p-value of 1 when it is not positive; pass it to do the
        same.  ``None`` ignores the rule.

    Raises
    ------
    SpecificationError
        If ``offsets`` does not describe the residuals.
    NumericalError
        If no unit has more than one observation or all residuals are zero.
    """
    # Stata: xttest0 -> LM = (N^2 / 2) * A1^2 / (sum(T_i^2) - N), A1 = 1 - sum_i(sum_t e)^2 /
    # sum(e^2); p = chi2tail(1, LM) / 2; LM = 0 and p = 1 when e(sigma_u) <= 0
    e = as_array("residuals", residuals, ndim=1)
    bounds = np.asarray(offsets, dtype=np.int64)
    if (
        bounds.ndim != 1
        or bounds.size < 2
        or bounds[0] != 0
        or bounds[-1] != e.size
        or np.any(np.diff(bounds) < 1)
    ):
        raise SpecificationError(
            _i18n.message("diag.lm_random_effects.offsets", n=e.size),
            key="diag.lm_random_effects.offsets",
        )
    t = np.diff(bounds).astype(np.float64)
    n = float(e.size)
    m11 = float(np.sum(t**2))
    if m11 - n <= 0:
        raise NumericalError(
            _i18n.message("diag.lm_random_effects.single"), key="diag.lm_random_effects.single"
        )
    sum_sq = float(e @ e)
    if sum_sq == 0.0:
        raise NumericalError(
            _i18n.message("diag.lm_random_effects.zero_residuals"),
            key="diag.lm_random_effects.zero_residuals",
        )
    unit_sums = np.add.reduceat(e, bounds[:-1])
    a1 = 1.0 - float(np.sum(unit_sums**2)) / sum_sq
    notes: tuple[Note, ...] = ()
    if variance_u is not None and variance_u <= 0:
        statistic, pvalue = 0.0, 1.0
        notes = (Note("diag.lm_random_effects.zero_effect"),)
    else:
        statistic = (n**2 / 2.0) * a1**2 / (m11 - n)
        pvalue = float(stats.chi2.sf(statistic, 1)) / 2.0
    return TestResult(
        name="lm_random_effects",
        statistic=statistic,
        pvalue=pvalue,
        statistic_name="chibar2(01)",
        df=(1.0,),
        n_obs=e.size,
        n_groups=int(t.size),
        notes=notes,
    )
