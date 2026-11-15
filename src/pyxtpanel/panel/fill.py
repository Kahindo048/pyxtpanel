"""Filling the calendar (Stata's ``tsfill``) and balancing a panel."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._types import BoolArray

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "fill",
    {
        "fill.method_invalid": {
            "en": "method must be 'drop' or 'largest_block', got {value!r}.",
            "fr": "method doit valoir 'drop' ou 'largest_block', reçu {value!r}.",
        },
        "fill.indicator_exists": {
            "en": "The indicator column '{name}' already exists in the data.",
            "fr": "La colonne indicatrice « {name} » existe déjà dans les données.",
        },
    },
)


def missing_pairs(panel: PanelBase, *, full: bool) -> tuple[np.ndarray, np.ndarray]:
    """Return the unit codes and ranks of the periods ``tsfill`` adds.

    Without ``full``, a unit gets the periods between its first and its last observation that
    it lacks; with ``full``, every period from the first to the last of the whole panel.
    """
    index = panel.index
    lo, hi = int(np.min(index.rank)), int(np.max(index.rank))
    first = np.full(index.n_units, lo) if full else index.rank[index.offsets[:-1]]
    last = np.full(index.n_units, hi) if full else index.rank[index.offsets[1:] - 1]
    length = last - first + 1
    unit = np.repeat(np.arange(index.n_units, dtype=np.int64), length)
    start = np.cumsum(length) - length
    rank = (
        np.arange(unit.size, dtype=np.int64) - np.repeat(start, length) + np.repeat(first, length)
    )
    span = hi - lo + 1
    present = np.zeros(index.n_units * span, dtype=np.bool_)
    present[index.unit * span + (index.rank - lo)] = True
    absent = ~present[unit * span + (rank - lo)]
    return unit[absent], rank[absent]


# Stata: tsfill, tsfill full
def tsfill(panel: PanelBase, *, full: bool, indicator: str | None) -> pd.DataFrame:
    """Return the data of ``panel`` with one added row per missing period (``tsfill``).

    The result is sorted by unit then period.  Added rows have the identifier and time filled and
    every other column missing.  Key columns that are index levels become columns.
    """
    base = panel.to_frame(order="sorted")
    keys = [n for n in (panel.id_name, panel.time_name) if n in base.index.names]
    if keys:
        base = base.reset_index(level=keys)
    if indicator is not None and indicator in base.columns:
        raise SpecificationError(
            _i18n.message("fill.indicator_exists", name=indicator),
            stata_code=110,
            key="fill.indicator_exists",
        )
    unit, rank = missing_pairs(panel, full=full)
    if unit.size:
        added = pd.DataFrame(
            {
                panel.id_name: panel.index.unit_labels[unit],
                panel.time_name: panel.axis.values(rank),
            }
        )
        added = added.reindex(columns=base.columns)
        base = pd.concat([base, added], ignore_index=True)
        flag = np.r_[np.zeros(len(base) - unit.size, dtype=bool), np.ones(unit.size, dtype=bool)]
    else:
        base = base.reset_index(drop=True)
        flag = np.zeros(len(base), dtype=bool)
    if indicator is not None:
        base[indicator] = flag
    return base.sort_values([panel.id_name, panel.time_name], kind="stable").reset_index(drop=True)


def balanced_mask(
    panel: PanelBase, method: Literal["drop", "largest_block"], min_obs: int | None
) -> BoolArray:
    """Mask over the sorted rows kept by :meth:`PanelBase.balance`.

    * ``"drop"`` keeps the units observed in every period of the panel.
    * ``"largest_block"`` keeps the window of consecutive periods and the units complete in it
      that together hold the most observations (ties: the longer window, then the earlier one).

    ``min_obs`` first removes the units with fewer observations.
    """
    if method not in ("drop", "largest_block"):
        raise SpecificationError(
            _i18n.message("fill.method_invalid", value=method),
            stata_code=198,
            key="fill.method_invalid",
        )
    index = panel.index
    unit_ok = np.ones(index.n_units, dtype=np.bool_)
    if min_obs is not None:
        unit_ok &= index.counts >= min_obs
    if method == "drop":
        unit_ok &= index.counts == np.unique(index.rank).size
        return unit_ok[index.unit]
    lo, hi = int(np.min(index.rank)), int(np.max(index.rank))
    span = hi - lo + 1
    present = np.zeros((index.n_units, span), dtype=np.bool_)
    present[index.unit, index.rank - lo] = True
    present &= unit_ok[:, None]
    best: tuple[int, int, int, int] = (-1, 0, 0, 0)  # (observations, length, -start, stop)
    window = (lo, lo)
    best_units = np.zeros(index.n_units, dtype=np.bool_)
    for start in range(span):
        alive = np.flatnonzero(present[:, start])
        for stop in range(start, span):
            alive = alive[present[alive, stop]]
            if alive.size == 0:
                break
            length = stop - start + 1
            score = (alive.size * length, length, -start, stop)
            if score > best:
                best, best_units = score, np.zeros(index.n_units, dtype=np.bool_)
                best_units[alive] = True
                window = (start + lo, stop + lo)
    low, high = window
    return best_units[index.unit] & (index.rank >= low) & (index.rank <= high)
