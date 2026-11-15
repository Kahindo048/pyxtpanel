"""Pesaran CD test of cross-sectional dependence on a panel variable (category A, L.3.1).

Validated against Stata 17: ``xtcsd, pesaran`` after ``xtreg, fe`` (residuals of a fixed-effects
fit on a strongly balanced panel and on a panel with entry, exit and holes, where every pair uses
its common periods) and ``xtcd2, noest`` on a variable of a balanced panel.  The cases are in
``validation/cases/diagnostics_crosssection_pesaran_cd.py``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import (
    Note,
    TestResult,
    TestSpec,
    panel_variables,
    register_test,
)
from pyxtpanel.diagnostics._shared.pairwise import pairwise_pearson

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.pesaran_cd.title": {
            "en": "Pesaran CD test of cross-sectional dependence",
            "fr": "Test CD de Pesaran de dépendance transversale",
        },
        "diag.pesaran_cd.null": {
            "en": "the cross-section correlations are zero on average (weak dependence)",
            "fr": "les corrélations entre unités sont nulles en moyenne (dépendance faible)",
        },
        "diag.pesaran_cd.alternative": {
            "en": "the units are cross-sectionally dependent",
            "fr": "les unités sont dépendantes entre elles",
        },
        "diag.pesaran_cd.pairwise": {
            "en": "each pair of units uses its common periods ({low} to {high}); the statistic "
            "weights the correlation of a pair by the square root of that number.",
            "fr": "chaque paire d'unités utilise ses périodes communes ({low} à {high}) ; la "
            "statistique pondère la corrélation d'une paire par la racine de ce nombre.",
        },
        "diag.pesaran_cd.min_overlap": {
            "en": "min_overlap must be at least 3, got {value!r}.",
            "fr": "min_overlap doit valoir au moins 3, reçu {value!r}.",
        },
    },
)


@register_test(
    "pesaran_cd",
    TestSpec(
        family="crosssection",
        category="A",
        wave=1,
        reference="Pesaran, M. H. (2015), Econometric Reviews 34(6-10), 1089-1117",
        stata="xtcsd; xtcd2",
        validation=("cd_balanced", "cd_unbalanced", "cd_gaps", "cd_variable", "cd_independent"),
        inputs="panel",
        aliases=("cd",),
        balanced="both",
        r="plm::pcdtest(test='cd')",
        summary="Test that the correlations between units are zero on average.",
    ),
)
def pesaran_cd(panel: PanelBase, variable: str, *, min_overlap: int = 3) -> TestResult:
    """Return the CD statistic ``sqrt(2 / (N (N - 1))) * sum_{i<j} sqrt(T_ij) * rho_ij``.

    ``rho_ij`` is the correlation of ``variable`` between units ``i`` and ``j`` over the ``T_ij``
    periods where both are observed.  Under the null of weak cross-sectional dependence the
    statistic is standard normal, whether the panel is balanced or not, and the p-value is
    two-sided.  Apply it to the residuals of a model by adding them to the frame before declaring
    the panel; the unit effects must have been removed first (fixed-effects residuals), as the
    statistic reads the correlation of the series around their means.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of a numeric variable; rows where it is missing are ignored.
    min_overlap
        Smallest number of common periods a pair of units may have (3 or more; Stata needs 3).

    Raises
    ------
    SpecificationError
        If ``min_overlap`` is below 3.
    DataError
        If a pair of units has fewer than ``min_overlap`` periods in common.
    NumericalError
        If a pair has no variation over its common periods.
    """
    # Stata: xtcsd, pesaran -> E_ij = corr(e_i, e_j) * sqrt(r(N)) over the common periods;
    # CD = sqrt(2 / (N (N - 1))) * (1' E 1 - trace E) / 2; p = 2 * (1 - normal(|CD|)).
    # xtcd2, noest gives the same CD on a balanced panel.
    # STATA-QUIRK: on an unbalanced panel xtcd2 does not use this formula (it fills the panel and
    # correlates the unit means with all non-missing values); its CD differs from xtcsd's, e.g.
    # 16.08 against 2.52 on the 'gaps' panel.  This function follows Pesaran (2015) and xtcsd.
    if min_overlap < 3:
        raise SpecificationError(
            _i18n.message("diag.pesaran_cd.min_overlap", value=min_overlap),
            key="diag.pesaran_cd.min_overlap",
        )
    arrays = panel_variables(panel, [variable], test="pesaran_cd", min_units=3)
    n = arrays.n_units
    correlations = pairwise_pearson(arrays, test="pesaran_cd", min_overlap=min_overlap)
    rho_pairs, overlap = correlations.pairs()
    pairs = rho_pairs.size
    statistic = math.sqrt(2.0 / (n * (n - 1.0))) * float(np.sum(np.sqrt(overlap) * rho_pairs))
    notes: tuple[Note, ...] = ()
    if np.min(overlap) != np.max(overlap):
        notes = (
            Note(
                "diag.pesaran_cd.pairwise",
                {"low": int(np.min(overlap)), "high": int(np.max(overlap))},
            ),
        )
    details = pd.DataFrame(
        {
            "value": [
                float(rho_pairs.mean()),
                float(np.abs(rho_pairs).mean()),
                float(pairs),
                float(np.min(overlap)),
                float(np.max(overlap)),
            ]
        },
        index=pd.Index(
            [
                "average_correlation",
                "average_abs_correlation",
                "pairs",
                "min_overlap",
                "max_overlap",
            ],
            name="quantity",
        ),
    )
    return TestResult(
        name="pesaran_cd",
        statistic=statistic,
        pvalue=2.0 * float(stats.norm.sf(abs(statistic))),
        statistic_name="z",
        df=(),
        n_obs=arrays.n_obs,
        n_groups=n,
        variable=variable,
        options={"min_overlap": min_overlap},
        details=details,
        notes=notes,
    )
