"""``plot(kind="attrition")``: units observed, entries and exits, and survival of units (F.2.8)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._types import Lang
from pyxtpanel.describe.attrition import attrition
from pyxtpanel.plot._data import time_values
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import (
    AQUA,
    BLUE,
    INK,
    MUTED,
    ORANGE,
    chart_style,
    concise_dates,
    finish,
    integer_ticks,
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
        "plot.attr.active": {"en": "Units observed", "fr": "Unités observées"},
        "plot.attr.entries": {"en": "entries", "fr": "entrées"},
        "plot.attr.exits": {"en": "exits", "fr": "sorties"},
        "plot.attr.flows": {"en": "Entries and exits", "fr": "Entrées et sorties"},
        "plot.attr.survival": {
            "en": "Share of units still observed",
            "fr": "Part des unités encore observées",
        },
        "plot.attr.life": {
            "en": "Periods since the first observation",
            "fr": "Périodes depuis la première observation",
        },
        "plot.attr.censored": {
            "en": "{n} units still observed in the last period (censored)",
            "fr": "{n} unités encore observées à la dernière période (censurées)",
        },
    },
)


@plot_kind(
    "attrition",
    PlotSpec(summary="Units observed by period with entries and exits, and Kaplan-Meier survival."),
)
def plot_attrition(panel: PanelBase, *, lang: Lang | None = None) -> Figure:
    """Draw the units observed in each period, the entries and exits, and the survival curve.

    The left panels use :meth:`PanelData.attrition`'s table by period; the right panel is its
    Kaplan-Meier curve of the observed life of the units, where a unit still observed in the last
    period is censored and not counted as an exit.
    """
    chosen = resolve_lang(lang)
    result = attrition(panel)
    table = result.by_period
    low = int(np.min(panel.index.rank))
    ranks = np.arange(low, low + len(table), dtype=np.int64)
    x = time_values(panel, ranks)
    with chart_style():
        fig, axes = new_figure(height=4.4, shape=(2, 2))
        top, bottom = axes[0, 0], axes[1, 0]
        side = fig.add_subplot(axes[0, 1].get_gridspec()[:, 1])
        axes[0, 1].remove()
        axes[1, 1].remove()

        top.plot(x, table["active"].to_numpy(), color=BLUE, linewidth=1.6)
        top.fill_between(x, 0, table["active"].to_numpy(), color=BLUE, alpha=0.12, linewidth=0)
        top.set_ylabel(label("plot.attr.active", chosen))
        top.set_ylim(bottom=0)
        top.tick_params(labelbottom=False)

        entries, exits = table["entries"].to_numpy(), -table["exits"].to_numpy()
        bottom.fill_between(x, 0, entries, step="mid", color=AQUA, alpha=0.85, linewidth=0)
        bottom.fill_between(x, 0, exits, step="mid", color=ORANGE, alpha=0.85, linewidth=0)
        bottom.axhline(0, color=MUTED, linewidth=0.6)
        bottom.set_ylabel(label("plot.attr.flows", chosen))
        for text, height, color, va in (
            ("plot.attr.entries", 0.93, AQUA, "top"),
            ("plot.attr.exits", 0.07, ORANGE, "bottom"),
        ):
            bottom.text(
                0.5,
                height,
                label(text, chosen),
                transform=bottom.transAxes,
                ha="center",
                va=va,
                fontsize=7.5,
                color=color,
            )
        bottom.set_xlabel(panel.time_name)
        bottom.sharex(top)
        if panel.axis.freq is not None:
            concise_dates(bottom)
        else:
            integer_ticks(bottom)

        curve = result.survival
        steps = np.concatenate([[0], curve.index.to_numpy()])
        survival = np.concatenate([[1.0], curve["survival"].to_numpy()])
        side.step(steps, survival, where="post", color=BLUE, linewidth=1.6)
        side.set_ylim(0, 1.02)
        side.set_xlabel(label("plot.attr.life", chosen))
        side.set_ylabel(label("plot.attr.survival", chosen))
        integer_ticks(side)
        n_censored = int(curve["censored"].sum())
        side.text(
            0.03,
            0.04,
            label("plot.attr.censored", chosen, n=n_censored),
            transform=side.transAxes,
            ha="left",
            va="bottom",
            fontsize=7.5,
            color=INK,
        )
        return finish(fig, chosen)
