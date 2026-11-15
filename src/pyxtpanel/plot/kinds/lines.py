"""``plot(kind="lines")``: trajectories of the units, like ``xtline`` (F.2.1)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._types import Lang
from pyxtpanel.ops.frame import term_columns
from pyxtpanel.plot._data import period_mean, sample_units, time_values, unit_lines
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import (
    AQUA,
    BLUE,
    GREY,
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
        "plot.lines.mean": {"en": "mean of all units", "fr": "moyenne de toutes les unités"},
        "plot.lines.sample": {
            "en": "{shown} of {total} units",
            "fr": "{shown} unités sur {total}",
        },
    },
)

_HIGHLIGHT = (BLUE, ORANGE, AQUA)


@plot_kind(
    "lines",
    PlotSpec(summary="Trajectories of a sample of units, with the mean by period.", stata="xtline"),
)
def plot_lines(
    panel: PanelBase,
    variable: str,
    *,
    n_units: int = 12,
    units: Sequence[Any] | None = None,
    seed: int = 0,
    highlight: Sequence[Any] = (),
    mean: bool = True,
    panels: bool = False,
    lang: Lang | None = None,
) -> Figure:
    """Draw the trajectory of each unit of a reproducible sample.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        Column or Stata-style term (``"D.x"``).
    n_units
        Units drawn when the panel has more (a random sample fixed by ``seed``).
    units
        Labels of the units to draw instead of a sample.
    seed
        Seed of the sample.
    highlight
        Labels of up to three units drawn in colour and named at the end of their line; the
        others are grey.
    mean
        Add the mean by period of all the units of the panel (not only of those drawn).
    panels
        One small panel per unit instead of superimposed lines.
    lang
        ``"en"`` or ``"fr"`` (default: the global option).
    """
    chosen = resolve_lang(lang)
    names, values = term_columns(panel, variable)
    y = values[:, 0]
    codes = sample_units(panel, units, n_units, seed)
    lines = unit_lines(panel, y, codes)
    with chart_style():
        if panels:
            columns = min(4, max(1, len(lines)))
            rows = -(-len(lines) // columns)
            fig, grid = new_figure(
                height=1.5 * rows + 0.6,
                shape=(rows, columns),
                sharex=True,
                sharey=True,
                squeeze=False,
            )
            flat = grid.ravel()
            for ax, line in zip(flat, lines, strict=False):
                ax.plot(line.x, line.y, color=BLUE, linewidth=1.0)
                ax.set_title(str(line.label), fontsize=8, color=MUTED)
            for ax in flat[len(lines) :]:
                ax.set_visible(False)
            flat[0].set_ylabel(names[0])
            if panel.axis.freq is not None:
                concise_dates(flat[0], max_ticks=3)
            else:
                integer_ticks(flat[0])
            return finish(fig, chosen)
        fig, ax = new_figure(height=3.2)
        marked = {str(h): i for i, h in enumerate(list(highlight)[: len(_HIGHLIGHT)])}
        for line in lines:
            colour = _HIGHLIGHT[marked[str(line.label)]] if str(line.label) in marked else GREY
            ax.plot(
                line.x,
                line.y,
                color=colour,
                linewidth=1.3 if colour != GREY else 0.8,
                alpha=1.0 if colour != GREY else 0.55,
                zorder=3 if colour != GREY else 2,
            )
            if colour != GREY:
                ok = np.flatnonzero(~np.isnan(line.y))
                if ok.size:
                    ax.annotate(
                        str(line.label),
                        (line.x[ok[-1]], line.y[ok[-1]]),
                        xytext=(3, 0),
                        textcoords="offset points",
                        fontsize=7.5,
                        color=colour,
                        va="center",
                    )
        if mean:
            periods, average = period_mean(panel, y)
            ax.plot(
                time_values(panel, periods),
                average,
                color=INK,
                linewidth=1.4,
                linestyle="--",
                zorder=4,
            )
            ok = np.flatnonzero(~np.isnan(average))
            if ok.size:
                ax.annotate(
                    label("plot.lines.mean", chosen),
                    (time_values(panel, periods)[ok[0]], average[ok[0]]),
                    xytext=(3, 6),
                    textcoords="offset points",
                    fontsize=7.5,
                    color=INK,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5},
                    zorder=5,
                )
        ax.set_xlabel(panel.time_name)
        ax.set_ylabel(names[0])
        if panel.axis.freq is not None:
            concise_dates(ax)
        else:
            integer_ticks(ax)
        if len(lines) < panel.n_units:
            ax.text(
                0.99,
                0.97,
                label("plot.lines.sample", chosen, shown=len(lines), total=panel.n_units),
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7.5,
                color=MUTED,
            )
        return finish(fig, chosen)
