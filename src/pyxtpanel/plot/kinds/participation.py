"""``plot(kind="participation")``: units by periods, sorted by pattern (F.2.4)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._types import IntArray, Lang
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import (
    BLUE,
    MUTED,
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
        "plot.part.units": {"en": "Units", "fr": "Unités"},
        "plot.part.summary": {
            "en": "{patterns} patterns, {share:.0f} % of the cells observed",
            "fr": "{patterns} motifs, {share:.0f} % des cellules observées",
        },
    },
)


@dataclass(frozen=True)
class ParticipationData:
    """Presence matrix with units sorted by pattern.

    Attributes
    ----------
    present
        ``(n_units, n_periods)`` boolean, rows in the order drawn (top to bottom).
    ranks
        Calendar rank of each column (every period from the first to the last of the panel).
    order
        Unit code of each row.
    n_patterns
        Distinct presence patterns.
    """

    present: np.ndarray
    ranks: IntArray
    order: IntArray
    n_patterns: int


def participation_data(panel: PanelBase) -> ParticipationData:
    """Presence of each unit in each calendar period, units grouped by pattern.

    The most frequent pattern comes first; within equal frequencies, patterns are ordered by
    their first observed period then lexicographically, so that the order is reproducible.
    """
    index = panel.index
    low, high = int(np.min(index.rank)), int(np.max(index.rank))
    present = np.zeros((index.n_units, high - low + 1), dtype=bool)
    present[index.unit, index.rank - low] = True
    rows, inverse, counts = np.unique(present, axis=0, return_inverse=True, return_counts=True)
    first = np.argmax(rows, axis=1)
    key = np.lexsort((*(rows[:, j] for j in range(rows.shape[1] - 1, -1, -1)), first, -counts))
    rank_of_group = np.empty(rows.shape[0], dtype=np.int64)
    rank_of_group[key] = np.arange(rows.shape[0])
    order = np.argsort(rank_of_group[inverse.ravel()], kind="stable")
    return ParticipationData(
        present=present[order],
        ranks=np.arange(low, high + 1, dtype=np.int64),
        order=order.astype(np.int64),
        n_patterns=int(rows.shape[0]),
    )


@plot_kind(
    "participation",
    PlotSpec(summary="Map of the periods in which each unit is observed, sorted by pattern."),
)
def plot_participation(panel: PanelBase, *, lang: Lang | None = None) -> Figure:
    """Draw a cell for each (unit, period) observed, with units grouped by participation pattern."""
    from matplotlib.colors import ListedColormap

    chosen = resolve_lang(lang)
    data = participation_data(panel)
    with chart_style():
        fig, ax = new_figure(height=3.4)
        ax.imshow(
            data.present,
            aspect="auto",
            interpolation="nearest",
            cmap=ListedColormap(["#f4f3ef", BLUE]),
            vmin=0,
            vmax=1,
        )
        ticks = np.arange(0, data.ranks.size, max(1, -(-data.ranks.size // 8)))
        ax.set_xticks(ticks, panel.axis.labels(data.ranks[ticks]))
        ax.grid(False)
        ax.set_xlabel(panel.time_name)
        ax.set_ylabel(label("plot.part.units", chosen))
        share = 100.0 * data.present.mean()
        ax.text(
            0.99,
            -0.2,
            label("plot.part.summary", chosen, patterns=data.n_patterns, share=share),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            color=MUTED,
        )
        return finish(fig, chosen)
