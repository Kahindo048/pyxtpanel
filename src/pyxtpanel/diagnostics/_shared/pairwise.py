"""Correlations between the units of a panel, pair by pair (L.3).

The cross-section dependence tests and the correlation matrix share two computations:

* :func:`pairwise_pearson`: the Pearson correlation of every pair of units over the periods where
  both are observed, with the number of common periods ``T_ij`` (the unbalanced-panel form of
  Pesaran's CD, of the Breusch-Pagan LM and of Stata's ``xtcsd``);
* :func:`common_period_ranks`: the Spearman correlations over the periods where *all* units are
  observed (Friedman's and Frees' tests, as ``xtcsd, friedman | frees`` does).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, NumericalError
from pyxtpanel._types import FloatArray, IntArray
from pyxtpanel.diagnostics._core.inputs import plain_label
from pyxtpanel.panel.arrays import PanelArrays

_i18n.register(
    "diag",
    {
        "diag.pairwise.overlap": {
            "en": "{test}: units {first!r} and {second!r} have {have} periods in common, fewer "
            "than the {need} the test needs; the panel is too unbalanced for it.",
            "fr": "{test} : les unités {first!r} et {second!r} ont {have} périodes en commun, "
            "moins que les {need} exigées ; le panel est trop déséquilibré pour ce test.",
        },
        "diag.pairwise.constant": {
            "en": "{test}: units {first!r} and {second!r} have no variation over their common "
            "periods, so their correlation is not defined.",
            "fr": "{test} : les unités {first!r} et {second!r} ne varient pas sur leurs "
            "périodes communes, leur corrélation n'est pas définie.",
        },
        "diag.pairwise.no_common": {
            "en": "{test}: no period has all {units} units observed (or fewer than {need} do); "
            "the test uses the periods common to every unit.",
            "fr": "{test} : aucune période (ou moins de {need}) n'a les {units} unités "
            "observées ; le test utilise les périodes communes à toutes les unités.",
        },
    },
)


@dataclass(frozen=True)
class PairwiseCorrelations:
    """Pearson correlations between units, over each pair's common periods.

    Attributes
    ----------
    rho
        ``(N, N)`` symmetric matrix with ones on the diagonal; ``NaN`` for a pair that is too
        short or has no variation (only when the computation was not strict).
    overlap
        ``(N, N)`` number of common periods of each pair (the number of periods of the unit on
        the diagonal).
    labels
        Original unit identifiers, in order.
    """

    rho: FloatArray
    overlap: IntArray
    labels: np.ndarray

    def pairs(self) -> tuple[FloatArray, IntArray]:
        """Return the correlations and overlaps of the pairs ``i < j`` (flat arrays)."""
        upper = np.triu_indices(self.rho.shape[0], k=1)
        return self.rho[upper], self.overlap[upper]


def _wide(arrays: PanelArrays) -> tuple[FloatArray, FloatArray]:
    """Return the values and the observed mask as ``(periods, units)`` arrays."""
    periods, row_time = np.unique(arrays.time_index, return_inverse=True)
    values = np.zeros((periods.size, arrays.n_units))
    seen = np.zeros((periods.size, arrays.n_units))
    values[row_time, arrays.unit_of_row] = arrays.values[:, 0]
    seen[row_time, arrays.unit_of_row] = 1.0
    return values, seen


def pairwise_pearson(
    arrays: PanelArrays, *, test: str, min_overlap: int = 3, strict: bool = True
) -> PairwiseCorrelations:
    """Return the Pearson correlation of every pair of units over their common periods.

    Parameters
    ----------
    arrays
        One variable (column 0), missing rows already removed.
    test
        Name shown in the messages.
    min_overlap
        Smallest number of common periods a pair may have.
    strict
        Raise when a pair has fewer common periods or no variation (the tests do); otherwise
        put ``NaN`` for that pair (the descriptive matrix does).

    Raises
    ------
    DataError
        If ``strict`` and a pair has fewer than ``min_overlap`` common periods.
    NumericalError
        If ``strict`` and a pair has no variation over its common periods.
    """
    n = arrays.n_units
    values, seen = _wide(arrays)
    # Centre each unit on its own mean: the correlations do not change and the sums stay small.
    means = np.bincount(arrays.unit_of_row, weights=arrays.values[:, 0], minlength=n)
    values = (values - means / arrays.counts) * seen
    common = seen.T @ seen
    labels = arrays.unit_labels
    upper = np.triu_indices(n, k=1)
    if strict:
        short = np.flatnonzero(common[upper] < min_overlap)
        if short.size:
            i, j = int(upper[0][short[0]]), int(upper[1][short[0]])
            raise DataError(
                _i18n.message(
                    "diag.pairwise.overlap",
                    test=test,
                    first=plain_label(labels[i]),
                    second=plain_label(labels[j]),
                    have=int(common[i, j]),
                    need=min_overlap,
                ),
                key="diag.pairwise.overlap",
            )
    with np.errstate(divide="ignore", invalid="ignore"):
        sum_i = values.T @ seen  # [i, j] = sum of x_i over the periods where j is observed
        sum_ii = (values**2).T @ seen
        variance = sum_ii - sum_i**2 / common
        covariance = values.T @ values - sum_i * sum_i.T / common
        rho = covariance / np.sqrt(variance * variance.T)
    rho[common < min_overlap] = np.nan
    flat = np.flatnonzero(~np.isfinite(rho[upper]))
    if strict and flat.size:
        i, j = int(upper[0][flat[0]]), int(upper[1][flat[0]])
        raise NumericalError(
            _i18n.message(
                "diag.pairwise.constant",
                test=test,
                first=plain_label(labels[i]),
                second=plain_label(labels[j]),
            ),
            key="diag.pairwise.constant",
        )
    np.fill_diagonal(rho, 1.0)
    return PairwiseCorrelations(rho=rho, overlap=common.astype(np.int64), labels=labels)


def common_period_ranks(
    arrays: PanelArrays, *, test: str, min_periods: int
) -> tuple[FloatArray, int]:
    """Return the Spearman correlations over the periods where every unit is observed.

    Returns
    -------
    (rho, periods)
        The ``(N, N)`` Spearman correlation matrix (ranks with ties averaged) and the number of
        common periods.

    Raises
    ------
    DataError
        If fewer than ``min_periods`` periods have all units observed.
    """
    values, seen = _wide(arrays)
    complete = seen.all(axis=1)
    periods = int(complete.sum())
    if periods < min_periods:
        raise DataError(
            _i18n.message(
                "diag.pairwise.no_common", test=test, units=arrays.n_units, need=min_periods
            ),
            key="diag.pairwise.no_common",
        )
    ranks = stats.rankdata(values[complete], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.atleast_2d(np.corrcoef(ranks, rowvar=False))
    flat = np.flatnonzero(~np.isfinite(rho[np.triu_indices_from(rho, k=1)]))
    if flat.size:  # a unit that never changes has tied ranks: its correlations are 0 / 0
        i, j = (int(v[flat[0]]) for v in np.triu_indices_from(rho, k=1))
        raise NumericalError(
            _i18n.message(
                "diag.pairwise.constant",
                test=test,
                first=plain_label(arrays.unit_labels[i]),
                second=plain_label(arrays.unit_labels[j]),
            ),
            key="diag.pairwise.constant",
        )
    return rho, periods
