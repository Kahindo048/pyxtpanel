"""Data prepared for the figures: plain arrays that the figure draws and the tests compare."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from pyxtpanel._types import FloatArray, IntArray

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase


def time_values(panel: PanelBase, ranks: IntArray) -> np.ndarray:
    """Values to put on the time axis: timestamps for dated panels, integers otherwise."""
    if panel.axis.freq is not None:
        return np.asarray(panel.axis.datetimes(ranks).to_numpy())
    return np.asarray(panel.axis.values(ranks).to_numpy())


@dataclass(frozen=True)
class UnitLine:
    """One unit's line: ``x`` and ``y`` with a break (``NaN``) wherever periods are missing."""

    label: Any
    x: np.ndarray
    y: FloatArray


def unit_lines(panel: PanelBase, values: FloatArray, codes: Sequence[int]) -> list[UnitLine]:
    """Lines of the units ``codes``; a calendar gap becomes a break instead of a joined segment."""
    index = panel.index
    lines = []
    for code in codes:
        rows = slice(int(index.offsets[code]), int(index.offsets[code + 1]))
        rank, y = index.rank[rows], values[rows].copy()
        x = time_values(panel, rank)
        breaks = np.flatnonzero(np.diff(rank) > 1) + 1
        if breaks.size:
            x = np.insert(x, breaks, x[breaks - 1])
            y = np.insert(y, breaks, np.nan)
        lines.append(UnitLine(index.unit_labels[code], x, y))
    return lines


def sample_units(
    panel: PanelBase, units: Sequence[Any] | None, n_units: int, seed: int
) -> list[int]:
    """Return the unit codes to draw: the given labels, else a seeded sample of ``n_units``."""
    labels = panel.index.unit_labels
    if units is not None:
        return [int(c) for c in np.flatnonzero(np.isin(labels, list(units)))]
    if n_units >= panel.n_units:
        return list(range(panel.n_units))
    chosen = np.random.default_rng(seed).choice(panel.n_units, size=n_units, replace=False)
    return sorted(int(c) for c in chosen)


def period_mean(panel: PanelBase, values: FloatArray) -> tuple[IntArray, FloatArray]:
    """Mean of the observed values in each observed period (ranks and means)."""
    periods, code = np.unique(panel.index.rank, return_inverse=True)
    valid = ~np.isnan(values)
    total = np.bincount(code, weights=np.where(valid, values, 0.0), minlength=periods.size)
    count = np.bincount(code, weights=valid.astype(float), minlength=periods.size)
    with np.errstate(invalid="ignore", divide="ignore"):
        return periods, np.where(count > 0, total / count, np.nan)
