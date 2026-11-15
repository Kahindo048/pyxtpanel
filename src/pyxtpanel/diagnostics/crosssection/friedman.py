"""Friedman's test of cross-sectional independence (category A, L.3.3).

Validated against Stata 17: ``xtcsd, friedman`` after ``xtreg, fe`` on a balanced panel and on
a panel with entry and exit, where Stata keeps the periods common to every unit; the within
residuals are exported from Stata and given to the function.  The cases are in
``validation/cases/diagnostics_crosssection_friedman.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.pairwise import common_period_ranks

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.friedman.title": {
            "en": "Friedman test of cross-sectional independence",
            "fr": "Test de Friedman d'indépendance entre unités",
        },
        "diag.friedman.null": {
            "en": "the units are cross-sectionally independent",
            "fr": "les unités sont indépendantes entre elles",
        },
        "diag.friedman.alternative": {
            "en": "the average rank correlation between units is not zero",
            "fr": "la corrélation de rang moyenne entre unités n'est pas nulle",
        },
        "diag.friedman.common": {
            "en": "the ranks use the {periods} periods where every unit is observed, and the "
            "statistic uses the shortest unit length ({shortest}), as xtcsd does.",
            "fr": "les rangs utilisent les {periods} périodes où toutes les unités sont "
            "observées, et la statistique utilise la plus courte longueur d'unité ({shortest}), "
            "comme xtcsd.",
        },
    },
)


@register_test(
    "friedman",
    TestSpec(
        family="crosssection",
        category="A",
        wave=2,
        reference="Friedman, M. (1937), J. of the American Statistical Association 32(200), "
        "675-701; Frees, E. W. (1995), J. of Econometrics 69(2), 393-414",
        stata="xtcsd",
        validation=("fri_balanced", "fri_trimmed"),
        inputs="panel",
        aliases=("xtcsd_friedman",),
        balanced="both",
        r="plm::pcdtest(test='rho') (related)",
        summary="Test that the average Spearman correlation between units is zero.",
    ),
)
def friedman(panel: PanelBase, variable: str) -> TestResult:
    """Return Friedman's statistic ``(T - 1) ((N - 1) r + 1)``, chi-square with ``N - 1`` df.

    ``r`` is the average Spearman rank correlation between the units over the periods where
    every unit is observed; ``T`` is the length of the shortest unit.  The test is valid for
    ``T`` large or small, and does not need normal errors.  An unbalanced panel is accepted as
    in Stata's ``xtcsd``: only the periods common to all units give the ranks.  Residuals are
    tested by adding them to the frame before declaring the panel.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of a numeric variable; rows where it is missing are ignored.

    Raises
    ------
    DataError
        If fewer than 3 periods have every unit observed.
    """
    # Stata: xtcsd, friedman -> Spearman rho_ij on the periods where every unit is observed;
    # r = mean of the off-diagonal rho; fried = (tb - 1) ((N - 1) r + 1), tb = e(g_min);
    # p = 1 - chi2(N - 1, fried).
    arrays = panel_variables(panel, [variable], test="friedman", min_units=3)
    n = arrays.n_units
    shortest = int(np.min(arrays.counts))
    rho, periods = common_period_ranks(arrays, test="friedman", min_periods=3)
    mean_rho = float(np.mean(rho[np.triu_indices(n, k=1)]))
    statistic = (shortest - 1.0) * ((n - 1.0) * mean_rho + 1.0)
    notes = (
        (Note("diag.friedman.common", {"periods": periods, "shortest": shortest}),)
        if periods != shortest
        else ()
    )
    details = pd.DataFrame(
        {"value": [mean_rho, float(periods), float(shortest)]},
        index=pd.Index(
            ["average_rank_correlation", "common_periods", "shortest_unit"], name="quantity"
        ),
    )
    return TestResult(
        name="friedman",
        statistic=float(statistic),
        pvalue=float(stats.chi2.sf(statistic, n - 1.0)),
        statistic_name="chi2",
        df=(n - 1.0,),
        n_obs=arrays.n_obs,
        n_groups=n,
        variable=variable,
        details=details,
        notes=notes,
    )
