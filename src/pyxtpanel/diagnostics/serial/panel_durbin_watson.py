"""Panel Durbin-Watson and Baltagi-Wu LBI statistics for AR(1) errors (category A, L.8.10).

Validated against Stata 17: ``xtregar, fe lbi`` (the modified Bhargava et al. Durbin-Watson
statistic ``e(d1)`` and the Baltagi-Wu LBI statistic ``e(LBI)``) on a panel with AR(1) errors, on
independent errors, on panels with holes and with entry and exit, on a panel observed every two
periods and on a panel with a unit seen once.  The cases are in
``validation/cases/diagnostics_serial_panel_durbin_watson.py``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.ols import is_perfect_fit

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.panel_durbin_watson.title": {
            "en": "Panel Durbin-Watson and Baltagi-Wu LBI statistics",
            "fr": "Statistiques de Durbin-Watson de panel et LBI de Baltagi-Wu",
        },
        "diag.panel_durbin_watson.null": {
            "en": "the idiosyncratic errors are not autocorrelated (rho = 0)",
            "fr": "les erreurs idiosyncratiques ne sont pas autocorrélées (rho = 0)",
        },
        "diag.panel_durbin_watson.alternative": {
            "en": "the idiosyncratic errors follow an AR(1) process (rho != 0)",
            "fr": "les erreurs idiosyncratiques suivent un processus AR(1) (rho != 0)",
        },
        "diag.panel_durbin_watson.no_regressor": {
            "en": "panel_durbin_watson needs at least one regressor besides the dependent "
            "variable.",
            "fr": "panel_durbin_watson exige au moins un régresseur en plus de la variable "
            "dépendante.",
        },
        "diag.panel_durbin_watson.zero_residuals": {
            "en": "the within residuals are all zero: the statistics are not defined.",
            "fr": "les résidus intra sont tous nuls : les statistiques ne sont pas définies.",
        },
        "diag.panel_durbin_watson.no_pvalue": {
            "en": "no p-value: the null distribution depends on the regressors and the panel. "
            "Values near 2 are compatible with no autocorrelation, values well below 2 point to "
            "positive autocorrelation; compare with the bounds of Bhargava, Franzini and "
            "Narendranathan (1982) for the Durbin-Watson statistic and of Baltagi and Wu (1999) "
            "for the LBI statistic.",
            "fr": "pas de p-value : la loi nulle dépend des régresseurs et du panel. Des valeurs "
            "proches de 2 sont compatibles avec l'absence d'autocorrélation, des valeurs "
            "nettement inférieures à 2 indiquent une autocorrélation positive ; comparez aux "
            "bornes de Bhargava, Franzini et Narendranathan (1982) pour la statistique de "
            "Durbin-Watson et de Baltagi et Wu (1999) pour la statistique LBI.",
        },
    },
)


@register_test(
    "panel_durbin_watson",
    TestSpec(
        family="serial",
        category="A",
        wave=2,
        reference="Bhargava, A., Franzini, L. and Narendranathan, W. (1982), Review of Economic "
        "Studies 49(4), 533-549; Baltagi, B. H. and Wu, P. X. (1999), Econometric Theory 15(6), "
        "814-823",
        stata="xtregar",
        validation=(
            "pd_autocorrelated",
            "pd_independent",
            "pd_gaps",
            "pd_entry_exit",
            "pd_delta2",
            "pd_singleton",
        ),
        inputs="panel",
        balanced="both",
        r="plm::pbnftest",
        aliases=("xtregar", "lbi", "bhargava"),
        summary="Modified Durbin-Watson and Baltagi-Wu LBI statistics of the within residuals.",
    ),
)
def panel_durbin_watson(panel: PanelBase, dependent: str, regressors: Sequence[str]) -> TestResult:
    """Return the modified Durbin-Watson statistic of the within residuals, and the LBI.

    The model is fitted by within (fixed-effects) least squares.  With ``e`` its residuals, ``S``
    their sum of squares, and consecutive periods meaning one step of the calendar (a hole breaks
    the chain), the modified Durbin-Watson statistic of Bhargava, Franzini and Narendranathan is
    ``d1 = sum (e_it - e_i,t-1)^2 / S`` over the pairs of consecutive periods, where an
    observation that follows a hole enters with its own square.  Baltagi and Wu's locally best
    invariant statistic adds the squares of the residual before each hole, of the first and of the
    last residual of each unit: ``LBI = d1 + d2 + d3 + d4``.  For a balanced panel this is the
    Durbin-Watson statistic ``d1`` and ``LBI = d1 + (e_i1^2 + e_iT^2) / S`` summed over units.

    Both statistics are about 2 without autocorrelation and fall towards 0 with positive
    autocorrelation.  The result has no p-value (the null distribution depends on the data): the
    result carries a remark pointing to the published bounds.  ``statistic`` is ``d1``;
    ``details`` gives ``d1`` and the LBI.

    Parameters
    ----------
    panel
        The declared panel, balanced or not.
    dependent
        Name of the dependent variable.
    regressors
        Names of the regressors, without a constant (the unit effects are removed by the within
        transformation).

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If there is no regressor.
    NumericalError
        If the within residuals are all zero.
    """
    # Stata: xtregar, fe lbi -> lbi2/lbiw2: demean by id, _regress; Itj = (t - t[_n-1] == delta);
    # d1 = sum_{n>1} (e - e[_n-1] * Itj)^2 / S; d2 = sum_{n<N} e^2 (1 - Itjf) / S; d3, d4 = share
    # of e^2 at the first and the last observation of each unit; LBI = d1 + d2 + d3 + d4
    names = list(regressors)
    if not names:
        raise SpecificationError(
            _i18n.message("diag.panel_durbin_watson.no_regressor"),
            key="diag.panel_durbin_watson.no_regressor",
        )
    arrays = panel_variables(panel, [dependent, *names], test="panel_durbin_watson")
    unit = arrays.unit_of_row
    counts = np.bincount(unit)
    centred = (
        arrays.values
        - (np.add.reduceat(arrays.values, arrays.offsets[:-1], axis=0) / counts[:, None])[unit]
    )
    y, x = centred[:, 0], centred[:, 1:]
    coefficients, *_ = np.linalg.lstsq(x, y, rcond=1e-10)
    e = y - x @ coefficients
    total = float(e @ e)
    if is_perfect_fit(total, float(y @ y)):  # rounding noise is no residual
        raise NumericalError(
            _i18n.message("diag.panel_durbin_watson.zero_residuals"),
            key="diag.panel_durbin_watson.zero_residuals",
        )
    n = arrays.n_obs
    first = np.zeros(n, dtype=bool)
    first[arrays.offsets[:-1]] = True
    last = np.zeros(n, dtype=bool)
    last[arrays.offsets[1:] - 1] = True
    step = np.zeros(n, dtype=bool)  # the observation follows the previous one by one period
    step[1:] = (arrays.time_index[1:] - arrays.time_index[:-1] == 1) & ~first[1:]
    previous = np.zeros(n)
    previous[1:] = e[:-1]
    d1 = float(np.sum(((e - previous * step) ** 2)[~first])) / total
    followed = np.zeros(n, dtype=bool)  # the next observation is the next period
    followed[:-1] = step[1:]
    d2 = float(np.sum((e**2 * ~followed)[~last])) / total
    d3 = float(np.sum(e[first] ** 2)) / total
    d4 = float(np.sum(e[last] ** 2)) / total
    lbi = d1 + d2 + d3 + d4
    details = pd.DataFrame(
        {"value": [d1, lbi]}, index=pd.Index(["durbin_watson", "baltagi_wu_lbi"], name="statistic")
    )
    return TestResult(
        name="panel_durbin_watson",
        statistic=d1,
        pvalue=float("nan"),
        statistic_name="d",
        n_obs=n,
        n_groups=arrays.n_units,
        variable=dependent,
        options={"regressors": tuple(names)},
        details=details,
        notes=(Note("diag.panel_durbin_watson.no_pvalue"),),
    )
