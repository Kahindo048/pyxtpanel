"""Gaps inside units and repeated (unit, period) pairs in raw data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel.panel.declare import _levels

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase


def gaps_table(panel: PanelBase, *, only_with_gaps: bool = False) -> pd.DataFrame:
    """Per unit: first and last period, observations, breaks and missing periods inside.

    A *break* is an observation that follows the unit's previous one by more than one period.
    A *missing period* is a period between the unit's first and last one with no observation.

    Parameters
    ----------
    panel
        A declared panel.
    only_with_gaps
        Keep only the units that have at least one break.

    Returns
    -------
    pandas.DataFrame
        One row per unit, indexed by the unit identifier, with columns ``first``, ``last``,
        ``n_obs``, ``n_breaks`` and ``n_missing``.
    """
    index = panel.index
    starts, stops = index.offsets[:-1], index.offsets[1:]
    break_at = np.zeros(index.n_obs, dtype=np.int64)
    if index.n_obs > 1:
        same_unit = index.unit[1:] == index.unit[:-1]
        break_at[1:] = same_unit & (np.diff(index.rank) != 1)
    first, last = index.rank[starts], index.rank[stops - 1]
    table = pd.DataFrame(
        {
            "first": panel.axis.values(first),
            "last": panel.axis.values(last),
            "n_obs": index.counts,
            "n_breaks": np.add.reduceat(break_at, starts) if index.n_obs else break_at,
            "n_missing": (last - first + 1) - index.counts,
        },
        index=pd.Index(index.unit_labels, name=panel.id_name),
    )
    return table[table["n_breaks"] > 0] if only_with_gaps else table


def find_duplicates(
    frame: pd.DataFrame,
    id: str | None = None,
    time: str | None = None,  # noqa: A002
) -> pd.DataFrame:
    """Rows of ``frame`` whose (identifier, time) pair appears more than once.

    Use it on raw data that :func:`pyxtpanel.xtset` refused, to see what to merge or drop.

    Parameters
    ----------
    frame
        The raw data.
    id
        Column that identifies the unit (``None`` when a two-level index holds the keys).
    time
        Column that holds the period (``None`` when a two-level index holds the keys).

    Returns
    -------
    pandas.DataFrame
        The repeated rows with their original index, sorted by the two keys, and a column
        ``n_repeats`` giving how many rows share the pair.  Rows with a missing key are ignored.
    """
    id_series, time_series = _levels(frame, id, time)
    keys = pd.DataFrame({"i": id_series.to_numpy(), "t": time_series.to_numpy()})
    usable = keys.notna().all(axis=1).to_numpy()
    repeated = keys[keys.duplicated(keep=False).to_numpy() & usable]
    sizes = repeated.groupby(["i", "t"], sort=False)["i"].transform("size")
    ordered = repeated.sort_values(["i", "t"], kind="stable").index.to_numpy()
    out = frame.iloc[ordered].copy()
    out["n_repeats"] = sizes.loc[ordered].to_numpy()
    return out
