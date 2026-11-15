"""Typology of the gaps of an unbalanced panel (K.2.1).

An addition of the package (Stata, R and Python only detect gaps); validated by hand on small
panels and against an independent loop implementation.

The calendar of the panel runs from its first to its last period (every step of ``delta``).
Each unit is then described by what happens outside and inside its own span:

* ``late_entry``: first observed after the first period of the panel;
* ``early_exit``: last observed before the last period of the panel;
* ``internal_gap``: absent from at least one period between its first and last observation.

A unit can have several of these.  An internal gap and a *re-entry* are the same event seen from
its two ends: the gap table gives, for each gap, the last period before it and the period of
re-entry.  A unit with none of the three is ``complete``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.explore._spec import ToolSpec, tool

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

#: Flags in the order they are joined to build the type of a unit.
FLAGS = ("late_entry", "early_exit", "internal_gap")

_i18n.register(
    "explore",
    {
        "explore.gaps.title": {
            "en": "Gap typology: {n_units} units, {n_periods} periods ({first} to {last})",
            "fr": "Typologie des trous : {n_units} unités, {n_periods} périodes ({first} à {last})",
        },
        "explore.gaps.lengths": {
            "en": "{n} internal gaps; length in periods: mean {mean:.2f}, longest {longest}",
            "fr": "{n} trous internes ; longueur en périodes : moyenne {mean:.2f}, "
            "maximale {longest}",
        },
        "explore.gaps.none": {"en": "No internal gap.", "fr": "Aucun trou interne."},
    },
)


@dataclass(frozen=True, kw_only=True)
class GapTypes(Result):
    """Result of :func:`gap_types`.

    Attributes
    ----------
    by_unit
        One row per unit: ``first``, ``last``, ``n_obs``, ``lead_missing`` (periods of the
        panel before its first observation), ``trail_missing`` (after its last), ``n_gaps``,
        ``gap_periods`` (periods missing inside), the three boolean flags and ``type``
        (flags joined with ``"+"``, or ``"complete"``).
    gaps
        One row per internal gap: ``unit``, ``after`` (last period observed before the gap),
        ``resume`` (period of re-entry) and ``length`` (periods missing).
    counts
        One row per ``type``: ``units`` and ``share`` of the units, most frequent first.
    length_distribution
        Number of internal gaps (``gaps``) for each gap length (``length``).
    n_periods
        Periods of the calendar of the panel.
    first, last
        First and last period of the panel (as the user's time values).
    """

    by_unit: pd.DataFrame
    gaps: pd.DataFrame
    counts: pd.DataFrame
    length_distribution: pd.DataFrame
    n_periods: int
    first: object
    last: object

    def to_frame(self) -> pd.DataFrame:
        """Return the number of units of each type."""
        return self.counts

    def summary(self, lang: Lang | None = None) -> str:
        """Return the counts by type and the distribution of the gap lengths as text."""
        n_units = len(self.by_unit)
        lines = [
            _i18n.message(
                "explore.gaps.title",
                lang,
                n_units=n_units,
                n_periods=self.n_periods,
                first=self.first,
                last=self.last,
            ),
            "",
            self.counts.to_string(float_format=lambda v: f"{v:.3f}"),
            "",
        ]
        if len(self.gaps):
            lengths = self.gaps["length"]
            lines.append(
                _i18n.message(
                    "explore.gaps.lengths",
                    lang,
                    n=len(self.gaps),
                    mean=float(lengths.mean()),
                    longest=int(lengths.max()),
                )
            )
            lines.append(self.length_distribution.to_string())
        else:
            lines.append(_i18n.message("explore.gaps.none", lang))
        return "\n".join(lines)


@tool(
    "gap_types",
    ToolSpec(
        summary="Late entry, early exit, internal gap and re-entry of each unit, gap lengths.",
        validation="hand-computed panels, independent loop implementation, number of gaps "
        "against Stata tsreport",
        stata="tsreport, panel detail (lists the gaps, does not classify them)",
        similar="Stata xtpattern (SSC), moderndid.scan_gaps",
    ),
)
def gap_types(panel: PanelBase) -> GapTypes:
    """Classify each unit's absences relative to the calendar of the panel.

    Parameters
    ----------
    panel
        A declared panel.

    Returns
    -------
    GapTypes
        The table by unit, the table of internal gaps and the counts by type.
    """
    index = panel.index
    low, high = int(np.min(index.rank)), int(np.max(index.rank))
    starts, stops = index.offsets[:-1], index.offsets[1:] - 1
    first, last = index.rank[starts], index.rank[stops]

    jump = np.zeros(index.n_obs, dtype=np.int64)
    if index.n_obs > 1:
        same_unit = index.unit[1:] == index.unit[:-1]
        jump[1:] = np.where(same_unit, np.diff(index.rank) - 1, 0)
    at = np.flatnonzero(jump > 0)
    gaps = pd.DataFrame(
        {
            "unit": index.unit_labels[index.unit[at]],
            "after": panel.axis.values(index.rank[at - 1]),
            "resume": panel.axis.values(index.rank[at]),
            "length": jump[at],
        }
    )

    n_gaps = np.bincount(index.unit[at], minlength=index.n_units)
    flags = {
        "late_entry": first > low,
        "early_exit": last < high,
        "internal_gap": n_gaps > 0,
    }
    kind = np.array(
        [
            "+".join(name for name in FLAGS if flags[name][i]) or "complete"
            for i in range(index.n_units)
        ],
        dtype=object,
    )
    by_unit = pd.DataFrame(
        {
            "first": panel.axis.values(first),
            "last": panel.axis.values(last),
            "n_obs": index.counts,
            "lead_missing": first - low,
            "trail_missing": high - last,
            "n_gaps": n_gaps,
            "gap_periods": np.bincount(index.unit[at], weights=jump[at], minlength=index.n_units)
            .round()
            .astype(np.int64),
            **flags,
            "type": kind,
        },
        index=pd.Index(index.unit_labels, name=panel.id_name),
    )
    tally = by_unit["type"].value_counts(sort=False)
    tally = tally.sort_index().sort_values(ascending=False, kind="stable")
    counts = pd.DataFrame(
        {"units": tally.to_numpy(), "share": tally.to_numpy() / index.n_units},
        index=pd.Index(tally.index, name="type"),
    )
    lengths = gaps["length"].value_counts().sort_index()
    return GapTypes(
        by_unit=by_unit,
        gaps=gaps,
        counts=counts,
        length_distribution=pd.DataFrame(
            {"gaps": lengths.to_numpy()}, index=pd.Index(lengths.index, name="length")
        ),
        n_periods=high - low + 1,
        first=panel.axis.values(np.array([low]))[0],
        last=panel.axis.values(np.array([high]))[0],
    )
