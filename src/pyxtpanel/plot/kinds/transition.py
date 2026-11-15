"""``plot(kind="transition")``: heat map of the transition matrix of a variable (F.2.5)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal

from pyxtpanel import _i18n
from pyxtpanel._types import Lang
from pyxtpanel.describe.trans import xttrans
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import (
    BLUE,
    INK,
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
        "plot.trans.from": {"en": "Value at the earlier observation", "fr": "Valeur avant"},
        "plot.trans.to": {"en": "Value at the next observation", "fr": "Valeur suivante"},
        "plot.trans.note": {
            "en": "{n} transitions, % of each row",
            "fr": "{n} transitions, % de chaque ligne",
        },
        "plot.trans.dropped": {
            "en": "{n} spanning a calendar gap left out",
            "fr": "{n} enjambant un trou du calendrier écartées",
        },
        "plot.trans.bridged": {
            "en": "{n} bridge a calendar gap",
            "fr": "{n} enjambent un trou du calendrier",
        },
    },
)


@plot_kind(
    "transition",
    PlotSpec(summary="Heat map of the transition matrix between consecutive observations."),
)
def plot_transition(
    panel: PanelBase,
    variable: str,
    *,
    gaps: Literal["bridge", "strict"] = "bridge",
    bins: int | Sequence[float] | None = None,
    lang: Lang | None = None,
) -> Figure:
    """Draw the row percentages of ``xttrans`` as a heat map with the value in each cell.

    The options are those of :func:`~pyxtpanel.describe.trans.xttrans`; the number of
    transitions that bridge a calendar gap is written on the figure when there are some.
    """
    from matplotlib.colors import LinearSegmentedColormap

    chosen = resolve_lang(lang)
    result = xttrans(panel, variable, gaps=gaps, bins=bins)
    matrix = result.percent.to_numpy()
    with chart_style():
        height = max(3.0, 0.45 * matrix.shape[0] + 1.4)
        fig, ax = new_figure(height=height)
        image = ax.imshow(
            matrix,
            cmap=LinearSegmentedColormap.from_list("pyxtpanel", ["#ffffff", BLUE]),
            vmin=0,
            vmax=100,
            aspect="auto",
        )
        ax.grid(False)
        ax.set_xticks(range(matrix.shape[1]), [f"{v:g}" for v in result.percent.columns])
        ax.set_yticks(range(matrix.shape[0]), [f"{v:g}" for v in result.percent.index])
        ax.set_xlabel(label("plot.trans.to", chosen))
        ax.set_ylabel(label("plot.trans.from", chosen))
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                dark = matrix[i, j] > 55
                ax.text(
                    j,
                    i,
                    f"{matrix[i, j]:.1f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if dark else INK,
                )
        note = label("plot.trans.note", chosen, n=result.n_trans)
        if result.n_bridged:
            key = "plot.trans.dropped" if gaps == "strict" else "plot.trans.bridged"
            note += "; " + label(key, chosen, n=result.n_bridged)
        ax.set_title(note, fontsize=8, color=MUTED, loc="left")
        fig.colorbar(image, ax=ax, shrink=0.8, pad=0.02, label="%")
        return finish(fig, chosen)
