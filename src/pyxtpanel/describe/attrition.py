"""Entries, exits and re-entries of units over the calendar, and survival of units (E.1.3).

An addition of the package (Stata has no equivalent); validated by hand on small panels.

For each calendar period between the first and the last of the panel:

* ``active``: units observed;
* ``entries``: units observed for the first time;
* ``exits``: units observed for the last time;
* ``reentries``: units observed again after having been absent;
* ``in_gap``: units absent although observed before and after.

``survival`` is a Kaplan–Meier curve of the observed life of the units (periods from the first to
the last observation), where a unit still observed in the last period of the panel is censored,
not counted as an exit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.describe._spec import DescribeSpec, description

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase


@dataclass(frozen=True, kw_only=True)
class Attrition(Result):
    """Result of :func:`attrition`.

    Attributes
    ----------
    by_period
        One row per calendar period: ``active``, ``entries``, ``exits``, ``reentries``,
        ``in_gap``.
    survival
        One row per life length ``d`` (periods): ``at_risk``, ``exits`` (uncensored),
        ``censored`` and the Kaplan–Meier ``survival``.
    units_with_gaps
        Units that are absent from at least one period between their first and last observation.
    """

    by_period: pd.DataFrame
    survival: pd.DataFrame
    units_with_gaps: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table by period."""
        return self.by_period

    def summary(self, lang: Lang | None = None) -> str:
        """Return the table by period as text."""
        return self.by_period.to_string() + f"\n\nunits with gaps: {self.units_with_gaps}"


@description(
    "attrition",
    DescribeSpec(
        stata=None,
        summary="Entries, exits and re-entries of units by period, and unit survival.",
    ),
)
def attrition(panel: PanelBase) -> Attrition:
    """Count entries, exits and re-entries by calendar period and estimate unit survival."""
    index = panel.index
    low, high = int(np.min(index.rank)), int(np.max(index.rank))
    size = high - low + 1
    start, stop = index.offsets[:-1], index.offsets[1:] - 1
    first, last = index.rank[start] - low, index.rank[stop] - low
    active = np.bincount(index.rank - low, minlength=size)
    entries = np.bincount(first, minlength=size)
    exits = np.bincount(last, minlength=size)
    same_unit = index.unit[1:] == index.unit[:-1]
    returned = same_unit & (np.diff(index.rank) > 1)
    reentries = np.bincount(index.rank[1:][returned] - low, minlength=size)
    started = np.cumsum(entries)
    ended_before = np.concatenate([[0], np.cumsum(exits)[:-1]])
    in_gap = (started - ended_before) - active
    periods = np.arange(low, high + 1)
    by_period = pd.DataFrame(
        {
            "active": active,
            "entries": entries,
            "exits": exits,
            "reentries": reentries,
            "in_gap": in_gap,
        },
        index=pd.Index(panel.axis.labels(periods), name=panel.time_name),
    )
    life = last - first + 1
    censored_unit = last == size - 1
    longest = int(np.max(life))
    lengths = np.arange(1, longest + 1)
    at_risk = np.array([(life >= d).sum() for d in lengths])
    ended = np.bincount(life[~censored_unit], minlength=longest + 1)[1:]
    censored = np.bincount(life[censored_unit], minlength=longest + 1)[1:]
    hazard = np.divide(ended, at_risk, out=np.zeros(longest), where=at_risk > 0)
    survival = pd.DataFrame(
        {
            "at_risk": at_risk,
            "exits": ended,
            "censored": censored,
            "survival": np.cumprod(1.0 - hazard),
        },
        index=pd.Index(lengths, name="periods"),
    )
    return Attrition(
        by_period=by_period,
        survival=survival,
        units_with_gaps=int(np.unique(index.unit[1:][returned]).size),
    )
