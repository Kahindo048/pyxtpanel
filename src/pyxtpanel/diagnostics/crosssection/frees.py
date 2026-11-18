"""Frees' test of cross-sectional independence (category A, L.3.3).

Validated against Stata 17: ``xtcsd, frees`` after ``xtreg, fe`` on a balanced panel (critical
values of the Q distribution), on a panel with entry and exit and on a long panel (normal
approximation and p-value); the within residuals are exported from Stata and given to the
function.  The cases are in ``validation/cases/diagnostics_crosssection_frees.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.pairwise import common_period_ranks
from pyxtpanel.diagnostics._shared.tables import load_table

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.frees.title": {
            "en": "Frees test of cross-sectional independence",
            "fr": "Test de Frees d'indépendance entre unités",
        },
        "diag.frees.null": {
            "en": "the units are cross-sectionally independent",
            "fr": "les unités sont indépendantes entre elles",
        },
        "diag.frees.alternative": {
            "en": "the average squared rank correlation between units exceeds 1 / (T - 1)",
            "fr": "la corrélation de rang moyenne au carré entre unités dépasse 1 / (T - 1)",
        },
        "diag.frees.critical": {
            "en": "no p-value for fewer than 31 common periods: compare Q with the critical "
            "values of Frees' Q distribution (10%: {cv10:.4f}, 5%: {cv5:.4f}, 1%: {cv1:.4f}).",
            "fr": "pas de p-value sous 31 périodes communes : comparez Q aux valeurs critiques de "
            "la loi Q de Frees (10 % : {cv10:.4f}, 5 % : {cv5:.4f}, 1 % : {cv1:.4f}).",
        },
        "diag.frees.common": {
            "en": "the ranks use the {periods} periods where every unit is observed, and the "
            "statistic uses the shortest unit length ({shortest}), as xtcsd does.",
            "fr": "les rangs utilisent les {periods} périodes où toutes les unités sont "
            "observées, et la statistique utilise la plus courte longueur d'unité ({shortest}), "
            "comme xtcsd.",
        },
    },
)


@register_test(
    "frees",
    TestSpec(
        family="crosssection",
        category="A",
        wave=2,
        reference="Frees, E. W. (1995), J. of Econometrics 69(2), 393-414",
        stata="xtcsd",
        validation=("frees_balanced", "frees_trimmed", "frees_long"),
        inputs="panel",
        tables=("frees_q.csv",),
        aliases=("xtcsd_frees",),
        balanced="both",
        r="plm::pcdtest(test='rho') (related)",
        summary="Test that the average squared Spearman correlation between units is 1 / (T - 1).",
    ),
)
def frees(panel: PanelBase, variable: str) -> TestResult:
    """Return Frees' statistic ``Q = N (r2 - 1 / (T - 1))`` with its critical values or p-value.

    ``r2`` is the average squared Spearman rank correlation between the units over the periods
    where every unit is observed, ``T`` the length of the shortest unit and ``m`` the number of
    common periods.  For ``m`` from 4 to 30 the result has no p-value: the critical values of
    Frees' Q distribution at 10%, 5% and 1% are in ``details`` (``cv_10``, ``cv_5``, ``cv_1``)
    and in a note.  For ``m >= 31`` the statistic is compared with a normal approximation.
    Residuals are tested by adding them to the frame before declaring the panel.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of a numeric variable; rows where it is missing are ignored.

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    DataError
        If fewer than 4 periods have every unit observed.
    """
    # Stata: xtcsd, frees -> rave2 = mean of the off-diagonal squared Spearman rho_ij on the
    # periods where every unit is observed; frees = N (rave2 - 1 / (tb - 1)), tb = e(g_min).
    # With m common periods, m < 31: critical values T4..T30 of Frees' Q (no p-value); else
    # seq = sqrt(32/25 (m+2)^2 / ((m-1)^3 (m+1)^2) + 4/25 (5m+6)^2 (m-3) / (m (m-1)^2 (m+1)^2))
    # and p = 2 (1 - normal(frees / seq)).
    # STATA-QUIRK: that p-value doubles the normal tail although Q is only large under the
    # alternative, and exceeds 1 for a negative Q; it is reproduced, capped at 1.
    arrays = panel_variables(panel, [variable], test="frees", min_units=3)
    n = arrays.n_units
    shortest = int(np.min(arrays.counts))
    rho, m = common_period_ranks(arrays, test="frees", min_periods=4)
    mean_square = float(np.mean(rho[np.triu_indices(n, k=1)] ** 2))
    statistic = n * (mean_square - 1.0 / (shortest - 1.0))
    notes: list[Note] = []
    cv = (math.nan, math.nan, math.nan)
    pvalue = math.nan
    if m < 31:
        table = load_table("frees_q.csv")
        row = table[table["t"] == m]
        values = dict(
            zip(row["level"].to_numpy(dtype=int), row["value"].to_numpy(dtype=float), strict=True)
        )
        cv = (float(values[1]), float(values[5]), float(values[10]))
        notes.append(Note("diag.frees.critical", {"cv10": cv[2], "cv5": cv[1], "cv1": cv[0]}))
    else:
        first = (32.0 / 25.0) * ((m + 2.0) ** 2 / ((m - 1.0) ** 3 * (m + 1.0) ** 2))
        second = (4.0 / 25.0) * (
            ((5.0 * m + 6.0) ** 2 * (m - 3.0)) / ((m * (m - 1.0) ** 2) * (m + 1.0) ** 2)
        )
        pvalue = min(1.0, 2.0 * float(stats.norm.sf(statistic / math.sqrt(first + second))))
    if m != shortest:
        notes.append(Note("diag.frees.common", {"periods": m, "shortest": shortest}))
    details = pd.DataFrame(
        {"value": [mean_square, float(m), float(shortest), cv[2], cv[1], cv[0]]},
        index=pd.Index(
            [
                "average_squared_rank_correlation",
                "common_periods",
                "shortest_unit",
                "cv_10",
                "cv_5",
                "cv_1",
            ],
            name="quantity",
        ),
    )
    return TestResult(
        name="frees",
        statistic=float(statistic),
        pvalue=pvalue,
        statistic_name="Q",
        df=(),
        n_obs=arrays.n_obs,
        n_groups=n,
        variable=variable,
        details=details,
        notes=tuple(notes),
    )
