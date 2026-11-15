"""Correlation matrix between the units of a panel (descriptive, L.3.4).

Not a test: it shows what the cross-section dependence tests summarise, the correlation of the
variable (or of residuals added to the frame) between every pair of units over the periods where
both are observed, with their average and average absolute value.  Validated against an
independent Mata computation in Stata 17 on a balanced panel, on a panel with entry, exit and
holes, and on a sparse panel where several pairs are too short; the cases are in
``validation/cases/diagnostics_crosssection_cross_correlation.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.diagnostics._core import panel_variables
from pyxtpanel.diagnostics._shared.pairwise import pairwise_pearson

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.cross_correlation.summary": {
            "en": "Correlation between units of {variable}: {units} units, {pairs} pairs; "
            "average {average:.4f}, average absolute value {absolute:.4f}; {undefined} pair(s) "
            "with fewer than {need} common periods or no variation are left out.",
            "fr": "Corrélation entre unités de {variable} : {units} unités, {pairs} paires ; "
            "moyenne {average:.4f}, valeur absolue moyenne {absolute:.4f} ; {undefined} "
            "paire(s) avec moins de {need} périodes communes ou sans variation écartée(s).",
        },
        "diag.cross_correlation.min_overlap": {
            "en": "min_overlap must be at least 2, got {value!r}.",
            "fr": "min_overlap doit valoir au moins 2, reçu {value!r}.",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class CrossCorrelation(Result):
    """Result of :func:`cross_correlation`.

    Attributes
    ----------
    variable
        The variable the correlations are computed on.
    correlation
        ``(N, N)`` correlation between units, indexed by unit; ``NaN`` for a pair with fewer than
        ``min_overlap`` common periods or without variation.
    overlap
        ``(N, N)`` number of common periods of each pair (the unit's own count on the diagonal).
    average, average_abs
        Mean of the correlations of the defined pairs and mean of their absolute values.
    n_pairs, n_undefined
        Number of pairs with a correlation and without.
    min_overlap
        The smallest number of common periods a pair needs.
    """

    variable: str
    correlation: pd.DataFrame
    overlap: pd.DataFrame
    average: float
    average_abs: float
    n_pairs: int
    n_undefined: int
    min_overlap: int

    def to_frame(self) -> pd.DataFrame:
        """Return the correlation matrix."""
        return self.correlation

    def summary(self, lang: Lang | None = None) -> str:
        """Return the averages as text (the matrix is ``result.correlation``)."""
        return _i18n.message(
            "diag.cross_correlation.summary",
            lang,
            variable=self.variable,
            units=len(self.correlation),
            pairs=self.n_pairs,
            average=self.average,
            absolute=self.average_abs,
            undefined=self.n_undefined,
            need=self.min_overlap,
        )


def cross_correlation(panel: PanelBase, variable: str, *, min_overlap: int = 3) -> CrossCorrelation:
    """Return the correlation of ``variable`` between every pair of units.

    Each pair uses the periods where both units are observed.  Unlike the tests, a pair that is
    too short or has no variation does not stop the computation: its correlation is ``NaN`` and
    it is left out of the averages (and counted in ``n_undefined``).

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Name of a numeric variable; rows where it is missing are ignored.  To look at the
        residuals of a model, add them to the frame before declaring the panel.
    min_overlap
        Smallest number of common periods for a defined correlation (2 or more).

    Raises
    ------
    SpecificationError
        If ``min_overlap`` is below 2.
    DataError
        If the panel has fewer than 2 units.
    """
    if min_overlap < 2:
        raise SpecificationError(
            _i18n.message("diag.cross_correlation.min_overlap", value=min_overlap),
            key="diag.cross_correlation.min_overlap",
        )
    arrays = panel_variables(panel, [variable], test="cross_correlation", min_units=2)
    found = pairwise_pearson(
        arrays, test="cross_correlation", min_overlap=min_overlap, strict=False
    )
    rho, _ = found.pairs()
    defined = np.isfinite(rho)
    labels = pd.Index([_plain(label) for label in found.labels], name="unit")
    return CrossCorrelation(
        variable=variable,
        correlation=pd.DataFrame(found.rho, index=labels, columns=labels),
        overlap=pd.DataFrame(found.overlap, index=labels, columns=labels),
        average=float(np.mean(rho[defined])) if defined.any() else math.nan,
        average_abs=float(np.mean(np.abs(rho[defined]))) if defined.any() else math.nan,
        n_pairs=int(defined.sum()),
        n_undefined=int((~defined).sum()),
        min_overlap=min_overlap,
    )


def _plain(label: object) -> object:
    return label.item() if isinstance(label, np.generic) else label
