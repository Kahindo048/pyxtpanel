"""Spells of a discrete state: durations and changes of each unit (K.3.2).

An addition of the package (the closest tool is the user-written Stata command ``tsspell``);
validated by hand on small panels, against an independent loop implementation and by invariance
tests.

A *spell* is a maximal run of consecutive calendar periods in which a unit has the same observed
value.  A spell ends, and the next one starts, at a change of value *or* at a break of the run:
a calendar gap or a missing value.  The two are told apart by ``start_reason`` and
``end_reason``:

* ``first`` / ``last``: the unit's first / last observation (the spell may continue outside what
  is observed, so its length is a lower bound: ``open_start`` / ``open_end``);
* ``change``: the neighbouring observation, in the adjacent period, has another value;
* ``gap``: the adjacent period is not observed or its value is missing.

A *change* of a unit is a spell that starts with ``change``; units with at least one are the
units that move between states.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.explore._spec import ToolSpec, tool
from pyxtpanel.panel.arrays import sorted_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "explore",
    {
        "explore.spells.title": {
            "en": "Spells of {variable}: {n} spells, {changing} of {units} units change state "
            "({changes} changes)",
            "fr": "Épisodes de {variable} : {n} épisodes, {changing} unités sur {units} changent "
            "d'état ({changes} changements)",
        },
        "explore.spells.by_state": {
            "en": "Length of the complete spells (neither end open), by state:",
            "fr": "Longueur des épisodes complets (aucune extrémité ouverte), par état :",
        },
        "explore.spells.too_many": {
            "en": "{name!r} has {count} distinct values (limit {limit}); spells need a discrete "
            "state, discretize it first or raise max_states.",
            "fr": "« {name} » a {count} valeurs distinctes (limite {limit}) ; les épisodes "
            "exigent un état discret, discrétisez-la ou augmentez max_states.",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class Spells(Result):
    """Result of :func:`spells`.

    Attributes
    ----------
    spells
        One row per spell: ``unit``, ``state``, ``start``, ``end``, ``length`` (periods),
        ``start_reason``, ``end_reason``, ``open_start`` and ``open_end``.
    by_unit
        One row per unit: ``n_obs``, ``n_spells``, ``n_changes`` and ``changes``.
    by_state
        One row per state: ``spells``, ``complete`` (neither end open), ``mean``, ``median``
        and ``max`` length of the complete spells.
    variable
        The variable examined.
    n_changes, n_units_changing, n_units
        Totals.
    """

    spells: pd.DataFrame
    by_unit: pd.DataFrame
    by_state: pd.DataFrame
    variable: str
    n_changes: int
    n_units_changing: int
    n_units: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table by state."""
        return self.by_state

    def summary(self, lang: Lang | None = None) -> str:
        """Return the totals and the lengths by state as text."""
        return "\n".join(
            [
                _i18n.message(
                    "explore.spells.title",
                    lang,
                    variable=self.variable,
                    n=len(self.spells),
                    changing=self.n_units_changing,
                    units=self.n_units,
                    changes=self.n_changes,
                ),
                "",
                _i18n.message("explore.spells.by_state", lang),
                self.by_state.to_string(float_format=lambda v: f"{v:.2f}"),
            ]
        )


@tool(
    "spells",
    ToolSpec(
        summary="Spells of a discrete state: duration, changes per unit, units that change.",
        validation="hand-computed panels, independent loop implementation, invariance tests",
        similar="Stata tsspell (SSC)",
    ),
)
def spells(panel: PanelBase, variable: str, *, max_states: int = 50) -> Spells:
    """Cut the history of each unit into spells of constant value of ``variable``.

    Parameters
    ----------
    panel
        A declared panel.
    variable
        A discrete variable (numeric, with at most ``max_states`` distinct values).
    max_states
        Refuse a variable with more distinct values, which is almost surely continuous.
    """
    index = panel.index
    x = sorted_columns(panel, [variable])[:, 0]
    valid = ~np.isnan(x)
    levels = np.unique(x[valid])
    if levels.size > max_states:
        raise SpecificationError(
            _i18n.message(
                "explore.spells.too_many", name=variable, count=levels.size, limit=max_states
            ),
            stata_code=134,
            key="explore.spells.too_many",
        )
    n = index.n_obs
    adjacent = np.zeros(n, dtype=np.bool_)  # row k follows row k-1 of the same unit, next period
    if n > 1:
        adjacent[1:] = (index.unit[1:] == index.unit[:-1]) & (np.diff(index.rank) == 1)
    same = np.zeros(n, dtype=np.bool_)
    if n > 1:
        same[1:] = valid[1:] & valid[:-1] & (x[1:] == x[:-1])
    continues = adjacent & same
    begins = valid & ~continues
    ends = np.zeros(n, dtype=np.bool_)
    ends[valid] = True
    ends[:-1] &= ~continues[1:] | ~valid[:-1]
    first_row = np.zeros(n, dtype=np.bool_)
    first_row[index.offsets[:-1]] = True
    last_row = np.zeros(n, dtype=np.bool_)
    last_row[index.offsets[1:] - 1] = True

    start, stop = np.flatnonzero(begins), np.flatnonzero(ends)
    previous = np.maximum(start - 1, 0)
    following = np.minimum(stop + 1, n - 1)
    start_reason = np.where(
        first_row[start],
        "first",
        np.where(adjacent[start] & valid[previous], "change", "gap"),
    )
    end_reason = np.where(
        last_row[stop],
        "last",
        np.where(adjacent[following] & valid[following], "change", "gap"),
    )
    state = x[start]
    if np.all(state == np.round(state)):
        state = state.astype(np.int64)
    table = pd.DataFrame(
        {
            "unit": index.unit_labels[index.unit[start]],
            "state": state,
            "start": panel.axis.values(index.rank[start]),
            "end": panel.axis.values(index.rank[stop]),
            "length": index.rank[stop] - index.rank[start] + 1,
            "start_reason": start_reason,
            "end_reason": end_reason,
            "open_start": start_reason == "first",
            "open_end": end_reason == "last",
        }
    )
    unit_of = index.unit[start]
    n_spells = np.bincount(unit_of, minlength=index.n_units)
    n_changes = np.bincount(unit_of[start_reason == "change"], minlength=index.n_units)
    by_unit = pd.DataFrame(
        {
            "n_obs": np.add.reduceat(valid.astype(np.int64), index.offsets[:-1]),
            "n_spells": n_spells,
            "n_changes": n_changes,
            "changes": n_changes > 0,
        },
        index=pd.Index(index.unit_labels, name=panel.id_name),
    )
    complete = table[~(table["open_start"] | table["open_end"])]
    grouped = complete.groupby("state")["length"]
    by_state = pd.DataFrame(
        {
            "spells": table.groupby("state")["length"].size(),
            "complete": grouped.size(),
            "mean": grouped.mean(),
            "median": grouped.median(),
            "max": grouped.max(),
        }
    ).sort_index()
    by_state["complete"] = by_state["complete"].fillna(0).astype(np.int64)
    by_state.index.name = "state"
    return Spells(
        spells=table,
        by_unit=by_unit,
        by_state=by_state,
        variable=variable,
        n_changes=int(n_changes.sum()),
        n_units_changing=int((n_changes > 0).sum()),
        n_units=index.n_units,
    )
