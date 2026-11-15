"""``plot(kind="balance_cost")``: what each way of balancing the panel keeps (F.2.8)."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pyxtpanel import _i18n
from pyxtpanel._types import Lang
from pyxtpanel.explore.balance import balance_cost
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import (
    BLUE,
    INK,
    MUTED,
    ORANGE,
    chart_style,
    direct_label,
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
        "plot.cost.units": {"en": "units", "fr": "unités"},
        "plot.cost.obs": {"en": "observations", "fr": "observations"},
        "plot.cost.share": {"en": "Share of the panel kept", "fr": "Part du panel conservée"},
        "plot.cost.k": {
            "en": "Minimum observations per unit, k",
            "fr": "Observations minimales par unité, k",
        },
        "plot.cost.length": {
            "en": "Length of the window (periods)",
            "fr": "Longueur de la fenêtre (périodes)",
        },
        "plot.cost.thresholds": {
            "en": "Keep units with T_i >= k",
            "fr": "Garder les unités avec T_i >= k",
        },
        "plot.cost.windows": {
            "en": "Balanced window of consecutive periods",
            "fr": "Fenêtre cylindrée de périodes consécutives",
        },
        "plot.cost.largest": {
            "en": "largest block: {units} units x {periods} periods",
            "fr": "plus grand bloc : {units} unités x {periods} périodes",
        },
    },
)


@plot_kind(
    "balance_cost",
    PlotSpec(summary="Share of units and observations kept by each minimum length and window."),
)
def plot_balance_cost(panel: PanelBase, *, lang: Lang | None = None) -> Figure:
    """Draw the cost of balancing: the minimum-length rule and the balanced windows.

    Left: units and observations kept when only units with at least ``k`` observations stay (the
    panel is still unbalanced).  Right: the best window of consecutive periods of each length,
    which gives a balanced panel; the largest balanced block (most observations) is marked.
    """
    chosen = resolve_lang(lang)
    cost = balance_cost(panel)
    with chart_style():
        fig, (left, right) = new_figure(height=3.3, shape=(1, 2), sharey=True)
        k = cost.thresholds.index.to_numpy()
        left.step(k, cost.thresholds["units_share"], where="post", color=BLUE, linewidth=1.6)
        left.step(k, cost.thresholds["obs_share"], where="post", color=ORANGE, linewidth=1.6)
        left.set_xlabel(label("plot.cost.k", chosen))
        left.set_ylabel(label("plot.cost.share", chosen))
        left.set_title(label("plot.cost.thresholds", chosen), fontsize=8, color=MUTED, loc="left")
        left.set_ylim(0, 1.03)
        integer_ticks(left)

        length = cost.windows.index.to_numpy()
        right.plot(length, cost.windows["units_share"], color=BLUE, linewidth=1.6)
        right.plot(length, cost.windows["obs_share"], color=ORANGE, linewidth=1.6)
        right.set_xlabel(label("plot.cost.length", chosen))
        right.set_title(label("plot.cost.windows", chosen), fontsize=8, color=MUTED, loc="left")
        integer_ticks(right)
        big = cost.largest
        right.scatter(
            [big.n_periods],
            [big.n_obs / cost.n_obs],
            color=INK,
            s=22,
            zorder=5,
            edgecolors="white",
            linewidths=0.8,
        )
        right.text(
            0.98,
            0.96,
            label("plot.cost.largest", chosen, units=big.n_units, periods=big.n_periods),
            transform=right.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            color=INK,
        )
        direct_label(left, float(k[-1]) * 0.55, 0.06, label("plot.cost.units", chosen), color=BLUE)
        direct_label(left, float(k[-1]) * 0.55, 0.0, label("plot.cost.obs", chosen), color=ORANGE)
        return finish(fig, chosen)
