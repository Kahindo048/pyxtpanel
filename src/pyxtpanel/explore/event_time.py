"""Time relative to an event, and how far each unit is observed around it (K.3.4).

An addition of the package (``pyfixest`` only ships simulated data for event studies); validated
by hand on small panels, against an independent loop implementation and by invariance tests.

The event of a unit is given either by a 0/1 *indicator* (the first period in which it is non-zero)
or by an *event period* column (the period itself, constant within the unit).  The relative time
of a row is its distance to the event in steps of the calendar (``delta``): 0 in the event
period, negative before, positive after.  Units without an event are *never treated*; their
relative time is missing.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import IntArray, Lang
from pyxtpanel.explore._spec import ToolSpec, tool
from pyxtpanel.panel.arrays import sorted_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "explore",
    {
        "explore.event.title": {
            "en": "Event time: {treated} of {units} units have an event ({never} never treated)",
            "fr": "Temps d'événement : {treated} unités sur {units} ont un événement "
            "({never} jamais traitées)",
        },
        "explore.event.bounds": {
            "en": "Every treated unit is observed at least {pre} periods before the event and "
            "{post} periods from it on (the earliest and the latest common bounds).",
            "fr": "Chaque unité traitée est observée au moins {pre} périodes avant l'événement "
            "et {post} périodes à partir de celui-ci (bornes communes).",
        },
        "explore.event.one": {
            "en": "give exactly one of indicator= and event_period=.",
            "fr": "donnez exactement l'un de indicator= et event_period=.",
        },
        "explore.event.varies": {
            "en": "event_period {name!r} takes several values within unit {unit!r}; it must be "
            "constant within a unit (missing for never treated).",
            "fr": "event_period « {name} » prend plusieurs valeurs dans l'unité {unit!r} ; elle "
            "doit être constante dans l'unité (manquante si jamais traitée).",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class EventTime(Result):
    """Result of :func:`event_time`.

    Attributes
    ----------
    relative
        Relative time of each row, on the index of the declared frame (missing for never
        treated units and for rows outside the panel).
    by_unit
        One row per unit: ``event`` (period of the event, missing if never treated), ``first``
        and ``last`` observed periods, ``min_rel`` and ``max_rel`` (relative time of the
        first and last observation), ``pre`` and ``post`` (observations before the event and
        from it on) and ``treated``.
    coverage
        One row per relative time ``k``: ``units`` observed at ``k`` and their ``share`` of the
        treated units.
    common_pre, common_post
        Periods before the event and from it on that every treated unit is observed within
        (the minimum over treated units of ``-min_rel`` and of ``max_rel``).
    n_treated, n_never, n_units
        Totals.
    """

    relative: pd.Series
    by_unit: pd.DataFrame
    coverage: pd.DataFrame
    common_pre: int
    common_post: int
    n_treated: int
    n_never: int
    n_units: int

    def to_frame(self) -> pd.DataFrame:
        """Return the coverage by relative time."""
        return self.coverage

    def summary(self, lang: Lang | None = None) -> str:
        """Return the totals, the common bounds and the coverage as text."""
        lines = [
            _i18n.message(
                "explore.event.title",
                lang,
                treated=self.n_treated,
                units=self.n_units,
                never=self.n_never,
            )
        ]
        if self.n_treated:
            lines.append(
                _i18n.message(
                    "explore.event.bounds", lang, pre=self.common_pre, post=self.common_post
                )
            )
            lines += ["", self.coverage.to_string(float_format=lambda v: f"{v:.3f}")]
        return "\n".join(lines)


@tool(
    "event_time",
    ToolSpec(
        summary="Time relative to an event for each row, and the periods each unit is observed "
        "around it.",
        validation="hand-computed panels, independent loop implementation, invariance tests",
    ),
)
def event_time(
    panel: PanelBase, *, indicator: str | None = None, event_period: str | None = None
) -> EventTime:
    """Compute the time relative to each unit's event.

    Parameters
    ----------
    panel
        A declared panel.
    indicator
        Name of a 0/1 variable; the event is the first period in which it is non-zero.
    event_period
        Name of a column holding the period of the event (a value of the time variable, constant
        within the unit, missing if the unit is never treated).

    Returns
    -------
    EventTime
        Result of :func:`event_time`.

    Raises
    ------
    SpecificationError
        If both or neither of ``indicator`` and ``event_period`` are given.
    DataError
        If ``event_period`` is not the same on all the rows of a unit.
    """
    if (indicator is None) == (event_period is None):
        raise SpecificationError(
            _i18n.message("explore.event.one"), stata_code=198, key="explore.event.one"
        )
    index = panel.index
    n_units = index.n_units
    never_rank = np.iinfo(np.int64).max
    event_rank = np.full(n_units, never_rank, dtype=np.int64)
    if indicator is not None:
        x = sorted_columns(panel, [indicator])[:, 0]
        on = (x != 0) & ~np.isnan(x)
        np.minimum.at(event_rank, index.unit[on], index.rank[on])
    else:
        assert event_period is not None
        column = panel.frame[event_period].iloc[index.source_row].reset_index(drop=True)
        given = column.notna().to_numpy()
        if given.any():
            ranks = np.zeros(index.n_obs, dtype=np.int64)
            ranks[given] = panel.axis.rank_of(column[given])
            low_u = np.full(n_units, never_rank, dtype=np.int64)
            high_u = np.full(n_units, -never_rank, dtype=np.int64)
            np.minimum.at(low_u, index.unit[given], ranks[given])
            np.maximum.at(high_u, index.unit[given], ranks[given])
            bad = np.flatnonzero((high_u != -never_rank) & (low_u != high_u))
            if bad.size:
                raise DataError(
                    _i18n.message(
                        "explore.event.varies",
                        name=event_period,
                        unit=index.unit_labels[bad[0]],
                    ),
                    stata_code=459,
                    key="explore.event.varies",
                )
            event_rank = np.where(high_u != -never_rank, low_u, never_rank)
    treated = event_rank != never_rank
    rel = np.where(treated[index.unit], index.rank - event_rank[index.unit], 0)
    sorted_rel = np.where(treated[index.unit], rel, np.nan)
    relative = cast("pd.Series", panel.align(sorted_rel, name="relative_time"))

    starts, stops = index.offsets[:-1], index.offsets[1:] - 1
    first, last = index.rank[starts], index.rank[stops]
    pre = np.add.reduceat((treated[index.unit] & (rel < 0)).astype(np.int64), starts)
    post = np.add.reduceat((treated[index.unit] & (rel >= 0)).astype(np.int64), starts)
    safe = np.where(treated, event_rank, 0)
    by_unit = pd.DataFrame(
        {
            "event": pd.Series(panel.axis.values(safe), index=range(n_units))
            .where(treated)
            .to_numpy(),
            "first": panel.axis.values(first),
            "last": panel.axis.values(last),
            "min_rel": np.where(treated, first - safe, np.nan),
            "max_rel": np.where(treated, last - safe, np.nan),
            "pre": np.where(treated, pre, 0),
            "post": np.where(treated, post, 0),
            "treated": treated,
        },
        index=pd.Index(index.unit_labels, name=panel.id_name),
    )
    n_treated = int(treated.sum())
    keep = treated[index.unit]
    steps: IntArray = rel[keep]
    if n_treated:
        ks, counts = np.unique(steps, return_counts=True)
        full = np.arange(int(np.min(ks)), int(np.max(ks)) + 1)
        units_at = np.zeros(full.size, dtype=np.int64)
        units_at[ks - int(np.min(ks))] = counts
        coverage = pd.DataFrame(
            {"units": units_at, "share": units_at / n_treated},
            index=pd.Index(full, name="k"),
        )
        common_pre = int(-by_unit.loc[treated, "min_rel"].max())
        common_post = int(by_unit.loc[treated, "max_rel"].min())
    else:
        coverage = pd.DataFrame({"units": [], "share": []}, index=pd.Index([], name="k"))
        common_pre = common_post = 0
    return EventTime(
        relative=relative,
        by_unit=by_unit,
        coverage=coverage,
        common_pre=max(common_pre, 0),
        common_post=max(common_post, 0),
        n_treated=n_treated,
        n_never=n_units - n_treated,
        n_units=n_units,
    )
