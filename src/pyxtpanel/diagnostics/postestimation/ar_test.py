"""Arellano-Bond test of serial correlation in the residuals, from bare arrays (L.6.1).

Validated against Stata 17: ``xtabond2`` (``e(ar1)``, ``e(ar2)`` and their p-values, one-step
and two-step robust) and ``xtabond`` followed by ``estat abond``; the cases are in
``validation/cases/diagnostics_postestimation_ar_test.py``.  Only the robust (clustered) form is
supported; the non-robust one-step form and the levels residuals of system GMM are refused by
the estimators' own code until they are validated.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import (
    Equations,
    Note,
    TestResult,
    TestSpec,
    as_array,
    register_test,
)

_i18n.register(
    "diag",
    {
        "diag.ar_test.title": {
            "en": "Arellano-Bond test for serial correlation in the residuals",
            "fr": "Test d'Arellano-Bond d'autocorrélation sérielle des résidus",
        },
        "diag.ar_test.null": {
            "en": "no autocorrelation of the given order in the differenced errors",
            "fr": "pas d'autocorrélation de l'ordre donné dans les erreurs en différence",
        },
        "diag.ar_test.alternative": {
            "en": "autocorrelation of the given order",
            "fr": "autocorrélation de l'ordre donné",
        },
        "diag.ar_test.no_pairs": {
            "en": "no pair of residuals {order} periods apart: the statistic is not defined.",
            "fr": "aucune paire de résidus distants de {order} périodes : la statistique "
            "n'est pas définie.",
        },
        "diag.ar_test.too_few": {
            "en": "{pairs} pairs from {units} units, fewer than the {minimum} required: the "
            "statistic is not computed.",
            "fr": "{pairs} paires issues de {units} unités, moins que les {minimum} requises : "
            "la statistique n'est pas calculée.",
        },
        "diag.ar_test.variance": {
            "en": "the variance of the numerator is not positive ({value:.6g}): the statistic "
            "is not defined.",
            "fr": "la variance du numérateur n'est pas positive ({value:.6g}) : la statistique "
            "n'est pas définie.",
        },
        "diag.ar_test.bad_argument": {
            "en": "{name} must be {expected}, got {value}.",
            "fr": "{name} doit être {expected}, reçu {value}.",
        },
        "diag.ar_test.mismatch": {
            "en": "the test and estimation systems must describe the same {what}: {a} against {b}.",
            "fr": "les systèmes de test et d'estimation doivent décrire les mêmes {what} : {a} "
            "contre {b}.",
        },
    },
)


@register_test(
    "ar_test",
    TestSpec(
        family="postestimation",
        category="B",
        wave=1,
        reference="Arellano, M. and Bond, S. (1991), Review of Economic Studies 58(2), 277-297",
        stata="xtabond2; estat abond",
        validation=(
            "ar_ab1_balanced",
            "ar_ab2_balanced",
            "ar_ab2_gaps",
            "ar_ab2_entry_exit",
            "ar_abond_balanced",
            "ar_abond1_balanced",
        ),
        inputs="arrays",
        balanced="na",
        r="plm::mtest",
        summary="Arellano-Bond m_k test of serial correlation of order k in differenced residuals.",
    ),
)
def ar_test(
    order: int,
    test: Equations,
    *,
    weight: object,
    covariance: object,
    estimation: Equations | None = None,
    min_pairs: int = 1,
) -> TestResult:
    """Return the Arellano-Bond statistic ``m_k`` of serial correlation of order ``order``.

    ``m_k = d0 / sqrt(v)`` with ``d0`` the sum over units of the products of the residuals of the
    test system with their value ``order`` periods earlier, and
    ``v = T1 + T2 + T3``: ``T1`` the sum over units of the squared unit sums of those products,
    ``T2 = -2 g' (X'Z W Z'X)^-1 X'Z W s`` the covariance with the instruments' moments and
    ``T3 = g' V g`` the effect of estimating the parameters, with
    ``g = sum_i w_i' X_i`` and ``s = sum_i d_i Z_i' e_i``.  The lag is read on the calendar, so a
    missing period gives no pair.

    Parameters
    ----------
    order
        Order ``k`` of the autocorrelation (1 or 2 in practice).
    test
        First-difference equations the test is run on (residuals, regressors, calendar).
    weight
        ``(L, L)`` weight matrix of the final estimation step.
    covariance
        ``(p, p)`` covariance matrix of the parameters (robust, Windmeijer-corrected after two
        steps).
    estimation
        Equations that were estimated, with their instruments, when they differ from ``test``
        (forward orthogonal deviations are estimated but tested in first differences); by default
        ``test``, whose instruments must then be given.
    min_pairs
        Fewest residual pairs, and units contributing one, for which the statistic is computed
        (``pyxtdpthresh`` uses 5).

    Returns
    -------
    TestResult
        ``statistic`` is ``m_k`` (a standard normal under the null), ``n_obs`` the number of
        pairs and ``n_groups`` the number of units that contribute; ``NaN`` with a remark when
        the statistic is not defined.

    Raises
    ------
    SpecificationError
        If an argument or the sizes of the systems do not match.
    """
    # STATA-QUIRK: on a panel with gaps xtabond and xtabond2 give different AR statistics (here
    # AR(1) -3.049 and -2.963); the equations handed to this function decide which one is reproduced
    # Stata: xtabond2 -> e(ar1), e(ar1p), e(ar2), e(ar2p); xtabond + estat abond -> r(arm), r(p)
    system = test if estimation is None else estimation
    _check_arguments(order, min_pairs, test, system)
    assert system.instruments is not None  # checked above
    z = system.instruments
    p = test.regressors.shape[1]
    w = as_array("weight", weight, ndim=2, shape=(z.shape[1], z.shape[1]))
    cov = as_array("covariance", covariance, ndim=2, shape=(p, p))

    unit = test.unit_of_row
    partner = _lag_partners(unit, test.time_index, order)
    found = partner >= 0
    lagged = np.where(found, test.residuals[np.maximum(partner, 0)], 0.0)
    unit_sum = np.bincount(unit, weights=test.residuals * lagged, minlength=test.n_units)
    n_pairs = int(found.sum())
    n_units = int(np.count_nonzero(np.bincount(unit[found], minlength=test.n_units)))
    d0, t1 = float(unit_sum.sum()), float(np.sum(unit_sum**2))

    notes: list[Note] = []
    statistic = math.nan
    t2 = t3 = math.nan
    if n_pairs == 0:
        notes.append(Note("diag.ar_test.no_pairs", {"order": order}))
    elif n_pairs < min_pairs or n_units < min_pairs:
        notes.append(
            Note(
                "diag.ar_test.too_few",
                {"pairs": n_pairs, "units": n_units, "minimum": min_pairs},
            )
        )
    else:
        g = test.regressors.T @ lagged
        s = z.T @ (system.residuals * unit_sum[system.unit_of_row])
        moments = system.regressors.T @ z  # (p, L)
        bread = np.linalg.pinv((moments @ w @ moments.T + (moments @ w @ moments.T).T) / 2.0)
        t2 = -2.0 * float(g @ (bread @ (moments @ (w @ s))))
        t3 = float(g @ cov @ g)
        variance = t1 + t2 + t3
        if variance > 0 and math.isfinite(variance):
            statistic = d0 / math.sqrt(variance)
        else:
            notes.append(Note("diag.ar_test.variance", {"value": variance}))
    pvalue = 2.0 * float(stats.norm.sf(abs(statistic))) if math.isfinite(statistic) else math.nan
    details = pd.DataFrame(
        {"value": [d0, t1, t2, t3]}, index=pd.Index(["d0", "T1", "T2", "T3"], name="term")
    )
    return TestResult(
        name="ar_test",
        statistic=statistic,
        pvalue=pvalue,
        statistic_name=f"m{order}",
        n_obs=n_pairs,
        n_groups=n_units,
        options={"order": order},
        details=details,
        notes=tuple(notes),
    )


def _check_arguments(order: int, min_pairs: int, test: Equations, system: Equations) -> None:
    for name, value, low in (("order", order, 1), ("min_pairs", min_pairs, 1)):
        if float(value) != int(value) or value < low:
            raise SpecificationError(
                _i18n.message(
                    "diag.ar_test.bad_argument",
                    name=name,
                    expected=f"an integer >= {low}",
                    value=value,
                ),
                key="diag.ar_test.bad_argument",
            )
    if system.instruments is None:
        raise SpecificationError(
            _i18n.message(
                "diag.ar_test.bad_argument",
                name="instruments",
                expected="given in the estimation system",
                value=None,
            ),
            key="diag.ar_test.bad_argument",
        )
    pairs = (
        ("units", test.n_units, system.n_units),
        ("regressors", test.regressors.shape[1], system.regressors.shape[1]),
    )
    for what, a, b in pairs:
        if a != b:
            raise SpecificationError(
                _i18n.message("diag.ar_test.mismatch", what=what, a=a, b=b),
                key="diag.ar_test.mismatch",
            )


def _lag_partners(unit: np.ndarray, time: np.ndarray, order: int) -> np.ndarray:
    """Row of the same unit ``order`` periods earlier, or -1 (rows sorted by unit and time)."""
    if time.size == 0:
        return np.empty(0, dtype=np.int64)
    low = int(np.min(time))
    span = int(np.max(time)) - low + order + 1
    key = unit.astype(np.int64) * span + (time - low + order)
    want = key - order
    at = np.minimum(np.searchsorted(key, want), key.size - 1)
    return np.where(key[at] == want, at, -1).astype(np.int64)
