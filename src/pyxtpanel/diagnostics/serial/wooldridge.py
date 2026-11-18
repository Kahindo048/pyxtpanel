"""Wooldridge test for first-order autocorrelation in panel data (category A, L.5.1).

Validated against Stata 17: ``xtserial`` (Drukker, Stata Journal 3(2)) on a panel with AR(1)
errors, on independent errors, and on panels with holes and with entry and exit, where the
differences and the lags follow the calendar.  The cases are in
``validation/cases/diagnostics_serial_wooldridge.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.ols import is_perfect_fit, ols_residuals

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.wooldridge.title": {
            "en": "Wooldridge test for autocorrelation in panel data",
            "fr": "Test de Wooldridge d'autocorrélation dans un panel",
        },
        "diag.wooldridge.null": {
            "en": "no first-order autocorrelation",
            "fr": "pas d'autocorrélation d'ordre un",
        },
        "diag.wooldridge.alternative": {
            "en": "the idiosyncratic errors are first-order autocorrelated",
            "fr": "les erreurs idiosyncratiques sont autocorrélées à l'ordre un",
        },
        "diag.wooldridge.no_regressor": {
            "en": "wooldridge needs at least one regressor besides the dependent variable.",
            "fr": "wooldridge exige au moins un régresseur en plus de la variable dépendante.",
        },
        "diag.wooldridge.few_pairs": {
            "en": "wooldridge: only {count} unit(s) have two consecutive first differences; the "
            "test needs at least 2 (three consecutive periods in a unit).",
            "fr": "wooldridge : seulement {count} unité(s) ont deux différences premières "
            "consécutives ; le test en exige au moins 2 (trois périodes consécutives dans une "
            "unité).",
        },
        "diag.wooldridge.no_variation": {
            "en": "wooldridge: the lagged first-difference residuals are all zero, so the "
            "autocorrelation is not defined.",
            "fr": "wooldridge : les résidus retardés en différences premières sont tous nuls, "
            "l'autocorrélation n'est pas définie.",
        },
        "diag.wooldridge.dropped": {
            "en": "{count} unit(s) lack two consecutive first differences and do not enter the "
            "second regression.",
            "fr": "{count} unité(s) n'ont pas deux différences premières consécutives et "
            "n'entrent pas dans la seconde régression.",
        },
    },
)


@register_test(
    "wooldridge",
    TestSpec(
        family="serial",
        category="A",
        wave=1,
        reference="Wooldridge, J. M. (2002), Econometric Analysis of Cross Section and Panel "
        "Data, MIT Press, sec. 10.4.1; Drukker, D. M. (2003), Stata Journal 3(2), 168-177",
        stata="xtserial",
        validation=("wr_autocorrelated", "wr_independent", "wr_gaps", "wr_entry_exit"),
        inputs="panel",
        aliases=("xtserial",),
        balanced="both",
        r="plm::pwfdtest",
        summary="Test of first-order serial correlation of the errors from first differences.",
    ),
)
def wooldridge(panel: PanelBase, dependent: str, regressors: Sequence[str]) -> TestResult:
    """Return the F test that the slope of the residuals on their lag is ``-1/2``.

    Under no serial correlation in the levels errors, the first-difference errors have a first
    autocorrelation of ``-1/2``.  The first differences of ``dependent`` are regressed on those of
    ``regressors`` without a constant; the residuals are then regressed, without a constant, on
    their own lag, and the slope is tested against ``-1/2`` with a covariance clustered by unit
    (the small-sample factor of Stata's ``regress, cluster``).  The statistic is
    ``F(1, G - 1)`` with ``G`` the number of units that enter the second regression.

    Differences and lags follow the calendar: a period missing inside a unit makes the
    differences on both sides of the hole missing, as in Stata, so the test accepts a panel with
    holes or with entry and exit.  Rows where any of the variables is missing are ignored.

    Parameters
    ----------
    panel
        The declared panel.
    dependent
        Name of the dependent variable of the model whose errors are tested.
    regressors
        Names of its regressors, without a constant (the first difference removes the unit effect
        and any constant).

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If there is no regressor.
    DataError
        If a variable is not numeric or no observation remains.
    NumericalError
        If fewer than two units have three consecutive periods, or the lagged residuals are zero.
    """
    # Stata: xtserial -> regress d.(y x), nocons cluster(id); predict res, residuals;
    # regress res l.res, nocons cluster(id); test l.res == -.5; F(1, G - 1), corr = _b[l.res].
    # The cluster variance is the small-sample one of regress: G / (G - 1) * (n - 1) / (n - 1).
    names = list(regressors)
    if not names:
        raise SpecificationError(
            _i18n.message("diag.wooldridge.no_regressor"), key="diag.wooldridge.no_regressor"
        )
    arrays = panel_variables(panel, [dependent, *names], test="wooldridge")
    previous = arrays.lag_index(1)
    has_diff = previous >= 0
    rows = np.flatnonzero(has_diff)
    delta = arrays.values[rows] - arrays.values[previous[rows]]
    residual = np.full(arrays.n_obs, np.nan)
    residual[rows] = ols_residuals(delta[:, 0], delta[:, 1:])
    # a perfect fit of the differences leaves rounding noise, not an exact zero
    if is_perfect_fit(float(residual[rows] @ residual[rows]), float(delta[:, 0] @ delta[:, 0])):
        raise NumericalError(
            _i18n.message("diag.wooldridge.no_variation"), key="diag.wooldridge.no_variation"
        )
    # Second regression: residual on its lag, for the rows where both exist (calendar lag).
    lagged = np.full(arrays.n_obs, np.nan)
    lagged[has_diff] = residual[previous[has_diff]]
    pair = np.flatnonzero(np.isfinite(residual) & np.isfinite(lagged))
    unit = arrays.unit_of_row[pair]
    clusters, cluster_code = np.unique(unit, return_inverse=True)
    groups = clusters.size
    if groups < 2:
        raise NumericalError(
            _i18n.message("diag.wooldridge.few_pairs", count=groups),
            key="diag.wooldridge.few_pairs",
        )
    x1, y1 = lagged[pair], residual[pair]
    sxx = float(x1 @ x1)
    if sxx == 0.0:
        raise NumericalError(
            _i18n.message("diag.wooldridge.no_variation"), key="diag.wooldridge.no_variation"
        )
    slope = float(x1 @ y1) / sxx
    score = np.bincount(cluster_code, weights=x1 * (y1 - slope * x1), minlength=groups)
    # Stata's regress, cluster(): (n - 1) / (n - k) * G / (G - 1) with k = 1 regressor.
    variance = groups / (groups - 1.0) * float(score @ score) / sxx**2
    statistic = (slope + 0.5) ** 2 / variance
    df_resid = float(groups - 1)
    left_out = arrays.n_units - groups
    notes = (Note("diag.wooldridge.dropped", {"count": left_out}),) if left_out else ()
    details = pd.DataFrame(
        {"value": [slope, math.sqrt(variance)]},
        index=pd.Index(["residual_autocorrelation", "standard_error"], name="quantity"),
    )
    return TestResult(
        name="wooldridge",
        statistic=statistic,
        pvalue=float(stats.f.sf(statistic, 1.0, df_resid)),
        statistic_name="F",
        df=(1.0, df_resid),
        n_obs=int(pair.size),
        n_groups=groups,
        variable=dependent,
        options={"regressors": tuple(names)},
        details=details,
        notes=notes,
    )
