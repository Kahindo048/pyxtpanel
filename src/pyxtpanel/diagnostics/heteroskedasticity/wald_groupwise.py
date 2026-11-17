"""Modified Wald test for groupwise heteroskedasticity, from bare arrays (category B, L.8.3).

Validated against Stata 17: ``xttest3`` (Baum, SSC) after ``xtreg, fe`` on a balanced panel with
unit-specific error variances and on an unbalanced panel; the within residuals and the unit
boundaries are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_heteroskedasticity_wald_groupwise.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, as_array, register_test

_i18n.register(
    "diag",
    {
        "diag.wald_groupwise.title": {
            "en": "Modified Wald test for groupwise heteroskedasticity",
            "fr": "Test de Wald modifié d'hétéroscédasticité par groupe",
        },
        "diag.wald_groupwise.null": {
            "en": "sigma(i)^2 = sigma^2 for all units i",
            "fr": "sigma(i)^2 = sigma^2 pour toutes les unités i",
        },
        "diag.wald_groupwise.alternative": {
            "en": "the error variance differs between units",
            "fr": "la variance des erreurs diffère selon les unités",
        },
        "diag.wald_groupwise.all_skipped": {
            "en": "no unit has a usable variance of its squared residuals (every unit has one "
            "observation or constant squared residuals): the test is not defined.",
            "fr": "aucune unité n'a de variance utilisable de ses résidus au carré (chaque unité a "
            "une observation ou des résidus au carré constants) : le test n'est pas défini.",
        },
        "diag.wald_groupwise.skipped": {
            "en": "{count} units have no usable variance of their squared residuals (a single "
            "period or an almost constant squared residual) and add nothing to the statistic, "
            "although they count in the degrees of freedom, as in Stata.",
            "fr": "{count} unités n'ont pas de variance exploitable de leurs résidus au carré "
            "(une seule période ou résidu au carré presque constant) et n'ajoutent rien à la "
            "statistique, bien qu'elles comptent dans les degrés de liberté, comme dans Stata.",
        },
        "diag.wald_groupwise.offsets": {
            "en": "offsets must start at 0, increase and end at the number of residuals ({n}).",
            "fr": "offsets doit commencer à 0, croître et finir au nombre de résidus ({n}).",
        },
    },
)

#: Smallest variance of the squared residuals of a unit that enters the statistic (as in xttest3).
MIN_VARIANCE = 1e-12


@register_test(
    "wald_groupwise",
    TestSpec(
        family="heteroskedasticity",
        category="B",
        wave=1,
        reference="Greene, W. H. (2000), Econometric Analysis, 4th ed., p. 598",
        stata="xttest3",
        validation=("wg_balanced", "wg_unbalanced", "wg_short_panel"),
        inputs="arrays",
        balanced="na",
        r="plm::pwtest (related); lmtest",
        summary="Test that the error variance is the same in every unit of a fixed-effects fit.",
    ),
)
def wald_groupwise(residuals: object, offsets: object) -> TestResult:
    """Return the modified Wald statistic ``sum_i (s_i^2 - s^2)^2 / V_i`` on ``G`` df.

    With ``e`` the within residuals, ``s_i^2`` the mean of ``e^2`` in unit ``i``, ``s^2`` the
    variance of all residuals and
    ``V_i = sum_t (e_it^2 - s_i^2)^2 / (T_i (T_i - 1))``, the statistic is a chi-square with as
    many degrees of freedom as units under the null of a common variance.  It is valid for
    unbalanced panels but **needs long panels**: ``V_i`` is an estimate of a fourth moment, so
    the test over-rejects when ``T_i`` is small (simulated with 40 units and normal errors, the
    nominal 5% test rejects always at ``T = 8``, about 37% of the time at ``T = 50`` and about
    10% at ``T = 200``).  Stata's ``xttest3`` behaves the same.

    Parameters
    ----------
    residuals
        ``(n,)`` residuals of the fixed-effects model (``e_it``, with the unit effect removed),
        sorted by unit.
    offsets
        ``(G + 1,)`` start of each unit in ``residuals`` (unit ``i`` owns
        ``residuals[offsets[i]:offsets[i + 1]]``).

    Raises
    ------
    SpecificationError
        If ``offsets`` does not describe the residuals.
    """
    # Stata: xttest3 -> s_i = mean_t(e^2); v_i = mean_t((e^2 - s_i)^2) / (T_i - 1), missing below
    # 1e-12; W = sum_i (s_i - var(e) * (N - 1) / N)^2 / v_i; chi2(e(N_g)) = W
    e = as_array("residuals", residuals, ndim=1)
    bounds = np.asarray(offsets, dtype=np.int64)
    if bounds.ndim != 1 or bounds.size < 2 or bounds[0] != 0 or bounds[-1] != e.size:
        raise SpecificationError(
            _i18n.message("diag.wald_groupwise.offsets", n=e.size),
            key="diag.wald_groupwise.offsets",
        )
    if np.any(np.diff(bounds) < 1):
        raise SpecificationError(
            _i18n.message("diag.wald_groupwise.offsets", n=e.size),
            key="diag.wald_groupwise.offsets",
        )
    groups = bounds.size - 1
    t = np.diff(bounds).astype(np.float64)
    unit = np.repeat(np.arange(groups), np.diff(bounds))
    e2 = e**2
    s_i = np.bincount(unit, weights=e2, minlength=groups) / t
    spread = np.bincount(unit, weights=(e2 - s_i[unit]) ** 2, minlength=groups) / t
    with np.errstate(divide="ignore", invalid="ignore"):
        v_i = np.where(t > 1, spread / (t - 1.0), np.nan)
    v_i = np.where(v_i < MIN_VARIANCE, np.nan, v_i)
    overall = float(np.mean((e - e.mean()) ** 2))
    contribution = (s_i - overall) ** 2 / v_i
    statistic = float(np.nansum(contribution))
    skipped = int(np.isnan(v_i).sum())
    if skipped == groups:
        raise NumericalError(
            _i18n.message("diag.wald_groupwise.all_skipped"),
            key="diag.wald_groupwise.all_skipped",
        )
    notes = (Note("diag.wald_groupwise.skipped", {"count": skipped}),) if skipped else ()
    details = pd.DataFrame(
        {"T": t.astype(int), "variance": s_i, "contribution": contribution},
        index=pd.RangeIndex(groups, name="unit"),
    )
    return TestResult(
        name="wald_groupwise",
        statistic=statistic if math.isfinite(statistic) else math.nan,
        pvalue=float(stats.chi2.sf(statistic, groups)),
        statistic_name="chi2",
        df=(float(groups),),
        n_obs=e.size,
        n_groups=groups,
        details=details,
        notes=notes,
    )
