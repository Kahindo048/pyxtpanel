"""Breusch-Pagan LM test of cross-sectional independence (category A, L.3.2).

Validated against Stata 17 with an independent Mata implementation of the published formula
(``xttest2``, the user-written command, is not installed): the within residuals of ``xtreg, fe``
are reshaped wide and the sum of ``T_ij * rho_ij^2`` over the pairs is computed in Mata, on a
panel with a common factor, on a panel with entry, exit and holes, and on independent series.
The cases are in ``validation/cases/diagnostics_crosssection_lm_cross_section.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core import Note, TestResult, TestSpec, panel_variables, register_test
from pyxtpanel.diagnostics._shared.pairwise import pairwise_pearson

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.lm_cross_section.title": {
            "en": "Breusch-Pagan LM test of cross-sectional independence",
            "fr": "Test LM de Breusch-Pagan d'indépendance entre unités",
        },
        "diag.lm_cross_section.null": {
            "en": "the units are cross-sectionally independent",
            "fr": "les unités sont indépendantes entre elles",
        },
        "diag.lm_cross_section.alternative": {
            "en": "some pairs of units are correlated",
            "fr": "certaines paires d'unités sont corrélées",
        },
        "diag.lm_cross_section.pairwise": {
            "en": "each pair of units uses its common periods ({low} to {high}); the statistic "
            "weights the squared correlation of a pair by that number.",
            "fr": "chaque paire d'unités utilise ses périodes communes ({low} à {high}) ; la "
            "statistique pondère le carré de la corrélation d'une paire par ce nombre.",
        },
        "diag.lm_cross_section.large_n": {
            "en": "the LM statistic is designed for few units relative to T; with N = {n} units "
            "and a shortest unit of {shortest} periods it over-rejects, use pesaran_cd.",
            "fr": "la statistique LM est conçue pour peu d'unités devant T ; avec N = {n} unités "
            "et une plus courte unité de {shortest} périodes elle rejette trop souvent, "
            "utilisez pesaran_cd.",
        },
        "diag.lm_cross_section.min_overlap": {
            "en": "min_overlap must be at least 3, got {value!r}.",
            "fr": "min_overlap doit valoir au moins 3, reçu {value!r}.",
        },
    },
)


@register_test(
    "lm_cross_section",
    TestSpec(
        family="crosssection",
        category="A",
        wave=2,
        reference="Breusch, T. S. and Pagan, A. R. (1980), Rev. Econ. Studies 47(1), 239-253",
        stata="mata",
        validation=("lmcs_factor", "lmcs_unbalanced", "lmcs_independent"),
        inputs="panel",
        aliases=("bp_lm", "xttest2"),
        balanced="both",
        r="plm::pcdtest(test='lm')",
        summary="Sum of the squared correlations between units, weighted by the common periods.",
    ),
)
def lm_cross_section(panel: PanelBase, variable: str, *, min_overlap: int = 3) -> TestResult:
    """Return ``LM = sum_{i<j} T_ij rho_ij^2``, chi-square with ``N (N - 1) / 2`` df.

    ``rho_ij`` is the correlation of ``variable`` between units ``i`` and ``j`` over their
    ``T_ij`` common periods.  The statistic suits panels with few units and many periods (the
    asymptotics are ``T -> infinity`` with ``N`` fixed); for a large ``N`` it over-rejects and
    ``pesaran_cd`` should be used.  Residuals are tested by adding them to the frame before
    declaring the panel, after the unit effects have been removed.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of a numeric variable; rows where it is missing are ignored.
    min_overlap
        Smallest number of common periods a pair of units may have (3 or more).

    Returns
    -------
    TestResult
        Outcome of a statistical test.

    Raises
    ------
    SpecificationError
        If ``min_overlap`` is below 3.
    DataError
        If a pair of units has fewer than ``min_overlap`` periods in common.
    NumericalError
        If a pair has no variation over its common periods.
    """
    # Stata: mata (no installed command; xttest2 of Baum computes the same): LM = sum_{i<j} T_ij
    # rho_ij^2 with rho_ij over the common periods of the wide residuals; chi2(N (N - 1) / 2).
    if min_overlap < 3:
        raise SpecificationError(
            _i18n.message("diag.lm_cross_section.min_overlap", value=min_overlap),
            key="diag.lm_cross_section.min_overlap",
        )
    arrays = panel_variables(panel, [variable], test="lm_cross_section", min_units=3)
    n = arrays.n_units
    correlations = pairwise_pearson(arrays, test="lm_cross_section", min_overlap=min_overlap)
    rho, overlap = correlations.pairs()
    statistic = float(np.sum(overlap * rho**2))
    df = n * (n - 1) / 2.0
    notes: list[Note] = []
    if np.min(overlap) != np.max(overlap):
        notes.append(
            Note(
                "diag.lm_cross_section.pairwise",
                {"low": int(np.min(overlap)), "high": int(np.max(overlap))},
            )
        )
    shortest = int(np.min(arrays.counts))
    if n > shortest:
        notes.append(Note("diag.lm_cross_section.large_n", {"n": n, "shortest": shortest}))
    details = pd.DataFrame(
        {"value": [float(np.mean(rho)), float(np.mean(np.abs(rho))), float(rho.size)]},
        index=pd.Index(
            ["average_correlation", "average_abs_correlation", "pairs"], name="quantity"
        ),
    )
    return TestResult(
        name="lm_cross_section",
        statistic=statistic,
        pvalue=float(stats.chi2.sf(statistic, df)),
        statistic_name="chi2",
        df=(df,),
        n_obs=arrays.n_obs,
        n_groups=n,
        variable=variable,
        options={"min_overlap": min_overlap},
        details=details,
        notes=tuple(notes),
    )
