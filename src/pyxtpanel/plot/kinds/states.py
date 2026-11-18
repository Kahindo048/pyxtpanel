"""``plot(kind="states")``: units by periods coloured by state, like ``panelView`` (F.2.7)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._types import IntArray, Lang
from pyxtpanel.panel.arrays import sorted_columns
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import (
    AQUA,
    BLUE,
    GREY,
    MUTED,
    ORANGE,
    chart_style,
    finish,
    label,
    new_figure,
    resolve_lang,
)

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "plot",
    {
        "plot.states.units": {"en": "Units", "fr": "Unités"},
        "plot.states.missing": {"en": "value missing", "fr": "valeur manquante"},
        "plot.states.absent": {"en": "not observed", "fr": "non observée"},
        "plot.states.changing": {
            "en": "{n} of {units} units change state",
            "fr": "{n} unités sur {units} changent d'état",
        },
        "plot.states.too_many": {
            "en": "{name!r} has {count} distinct values (limit {limit}); the map needs a "
            "discrete state, discretize it first or raise max_states.",
            "fr": "« {name} » a {count} valeurs distinctes (limite {limit}) ; la carte exige un "
            "état discret, discrétisez-la ou augmentez max_states.",
        },
        "plot.states.sort": {
            "en": "sort must be 'first_change' or 'unit', got {value!r}.",
            "fr": "sort doit valoir 'first_change' ou 'unit', reçu {value!r}.",
        },
    },
)

#: Colours of the states: a neutral first (the usual "untreated"), then the chart palette.
STATE_COLORS = ("#cdd5df", BLUE, ORANGE, AQUA, GREY, "#7a5bb8", "#c4a000", "#b3472d")
NOT_OBSERVED, MISSING_VALUE = "#ffffff", "#ebe9e3"

#: Codes of the cells of the grid below zero.
ABSENT, MISSING = -2, -1


@dataclass(frozen=True)
class StatesData:
    """The map as numbers.

    Attributes
    ----------
    grid
        ``(n_units, n_periods)`` integers, rows in the order drawn (top to bottom): ``-2`` the
        unit has no row in that period, ``-1`` it has a row but the value is missing, ``k`` the
        position of the state in ``states``.
    states
        Observed values of the variable, increasing.
    ranks
        Calendar rank of each column.
    order
        Unit code of each row.
    n_changing
        Units that take more than one state.
    """

    grid: np.ndarray
    states: np.ndarray
    ranks: IntArray
    order: IntArray
    n_changing: int


def states_data(
    panel: PanelBase,
    variable: str,
    *,
    sort: Literal["first_change", "unit"] = "first_change",
    max_states: int = 8,
) -> StatesData:
    """Grid of states of ``variable``, units sorted by the period of their first change.

    With ``sort="first_change"`` the units that change state earliest come first and the units
    that never change come last (in unit order); with ``"unit"`` the order is the unit order.
    """
    if sort not in ("first_change", "unit"):
        raise SpecificationError(
            _i18n.message("plot.states.sort", value=sort), stata_code=198, key="plot.states.sort"
        )
    index = panel.index
    x = sorted_columns(panel, [variable])[:, 0]
    valid = ~np.isnan(x)
    states = np.unique(x[valid])
    if states.size > max_states:
        raise SpecificationError(
            _i18n.message(
                "plot.states.too_many", name=variable, count=states.size, limit=max_states
            ),
            stata_code=134,
            key="plot.states.too_many",
        )
    low, high = int(np.min(index.rank)), int(np.max(index.rank))
    code = np.where(valid, np.searchsorted(states, np.where(valid, x, 0.0)), MISSING)
    grid = np.full((index.n_units, high - low + 1), ABSENT, dtype=np.int64)
    grid[index.unit, index.rank - low] = code

    rows = np.flatnonzero(valid)
    units, first = np.unique(index.unit[rows], return_index=True)
    first_code = np.full(index.n_units, -1, dtype=np.int64)
    first_code[units] = code[rows[first]]
    moved = valid & (code != first_code[index.unit])
    first_change = np.full(index.n_units, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(first_change, index.unit[moved], index.rank[moved])
    changing = first_change != np.iinfo(np.int64).max
    order = (
        np.lexsort((np.arange(index.n_units), first_change))
        if sort == "first_change"
        else np.arange(index.n_units)
    )
    return StatesData(
        grid=grid[order],
        states=states,
        ranks=np.arange(low, high + 1, dtype=np.int64),
        order=order.astype(np.int64),
        n_changing=int(changing.sum()),
    )


@plot_kind(
    "states",
    PlotSpec(
        summary="Map of units by periods coloured by the state of a discrete variable "
        "(treatment status), missing cells distinct."
    ),
)
def plot_states(
    panel: PanelBase,
    variable: str,
    *,
    sort: Literal["first_change", "unit"] = "first_change",
    max_states: int = 8,
    lang: Lang | None = None,
) -> Figure:
    """Draw one cell per (unit, period) coloured by the state of ``variable``.

    Cells where the unit has no row are white and cells where the value is missing are light
    grey, so that a missing treatment status is never mistaken for an untreated unit.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        A discrete numeric variable (a 0/1 treatment, a regime...).
    sort
        ``"first_change"`` (default) puts the units that change state earliest at the top, which
        shows the timing of a staggered treatment; ``"unit"`` keeps the order of the units.
    max_states
        Refuse a variable with more distinct values.
    lang
        ``"en"`` or ``"fr"`` (default: the global option).

    Returns
    -------
    matplotlib.figure.Figure
        The figure; save it with :func:`pyxtpanel.plot.save_figure`.
    """
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    chosen = resolve_lang(lang)
    data = states_data(panel, variable, sort=sort, max_states=max_states)
    k = data.states.size
    colors = [NOT_OBSERVED, MISSING_VALUE, *STATE_COLORS[:k]]
    with chart_style():
        fig, ax = new_figure(height=3.8)
        ax.imshow(
            data.grid + 2,
            aspect="auto",
            interpolation="nearest",
            cmap=ListedColormap(colors),
            vmin=0,
            vmax=k + 1,
        )
        ticks = np.arange(0, data.ranks.size, max(1, -(-data.ranks.size // 8)))
        ax.set_xticks(ticks, panel.axis.labels(data.ranks[ticks]))
        ax.grid(False)
        ax.set_xlabel(panel.time_name)
        ax.set_ylabel(label("plot.states.units", chosen))
        if panel.n_units <= 40:  # name the rows when they fit, otherwise leave them unlabelled
            ax.set_yticks(
                range(panel.n_units),
                [str(panel.index.unit_labels[c]) for c in data.order],
                fontsize=6.5,
            )
        else:
            ax.set_yticks([])
        handles = [
            Patch(facecolor=STATE_COLORS[i], edgecolor="none", label=f"{data.states[i]:g}")
            for i in range(k)
        ]
        if (data.grid == MISSING).any():
            handles.append(
                Patch(
                    facecolor=MISSING_VALUE,
                    edgecolor=MUTED,
                    label=label("plot.states.missing", chosen),
                )
            )
        if (data.grid == ABSENT).any():
            handles.append(
                Patch(
                    facecolor=NOT_OBSERVED,
                    edgecolor=MUTED,
                    label=label("plot.states.absent", chosen),
                )
            )
        ax.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.16),
            ncol=len(handles),
            fontsize=8,
            frameon=False,
        )
        ax.set_title(
            f"{variable}: "
            + label("plot.states.changing", chosen, n=data.n_changing, units=panel.n_units),
            fontsize=8,
            color=MUTED,
            loc="left",
        )
        return finish(fig, chosen)
