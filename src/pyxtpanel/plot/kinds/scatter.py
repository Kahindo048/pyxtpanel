"""``plot(kind="scatter")``: between and within scatter plots with a least-squares line (F.2.2)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._types import FloatArray, Lang
from pyxtpanel.describe.corr import unit_means
from pyxtpanel.ops.frame import term_columns
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import (
    BLUE,
    INK,
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
        "plot.scatter.between": {
            "en": "Between: unit means of {name}",
            "fr": "Entre : moyennes par unité de {name}",
        },
        "plot.scatter.within": {
            "en": "Within: deviations from the unit mean of {name}",
            "fr": "Intra : écarts à la moyenne de l'unité de {name}",
        },
        "plot.scatter.fit": {
            "en": "slope {slope:.3g}   n = {n}",
            "fr": "pente {slope:.3g}   n = {n}",
        },
        "plot.scatter.bad_part": {
            "en": "part must be 'between', 'within' or 'both', got {part!r}.",
            "fr": "part doit valoir 'between', 'within' ou 'both', reçu {part!r}.",
        },
    },
)


@dataclass(frozen=True)
class ScatterData:
    """Points and least-squares line of one panel of the figure."""

    x: FloatArray
    y: FloatArray
    slope: float
    intercept: float


def scatter_data(
    panel: PanelBase, x: FloatArray, y: FloatArray, kind: Literal["between", "within"]
) -> ScatterData:
    """Points of the between or within scatter over the rows where both variables are observed."""
    both = ~np.isnan(x) & ~np.isnan(y)
    unit = panel.index.unit[both]
    xs, ys = x[both], y[both]
    mx, present = unit_means(unit, panel.n_units, xs)
    my, _ = unit_means(unit, panel.n_units, ys)
    if kind == "between":
        px, py = mx[present], my[present]
    else:
        px, py = xs - mx[unit], ys - my[unit]
    if px.size > 1 and np.ptp(px) > 0:
        slope, intercept = (float(v) for v in np.polyfit(px, py, 1))
    else:
        slope, intercept = float("nan"), float("nan")
    return ScatterData(px, py, slope, intercept)


@plot_kind(
    "scatter",
    PlotSpec(summary="Between (unit means) and within (deviations) scatter plots of x and y."),
)
def plot_scatter(
    panel: PanelBase,
    x: str,
    y: str,
    *,
    part: Literal["between", "within", "both"] = "both",
    lang: Lang | None = None,
) -> Figure:
    """Scatter plot of ``y`` on ``x`` between units, within units, or both side by side.

    The between panel has one point per unit (the means over the rows where both variables are
    observed); the within panel has one point per observation (deviations from those means).
    Each panel carries its least-squares line and the slope and number of points, written on the
    figure instead of a legend.

    Parameters
    ----------
    panel
        The declared panel.
    x
        Name of the variable on the horizontal axis.
    y
        Name of the variable on the vertical axis.
    part
        ``"between"``, ``"within"`` or ``"both"`` (default): the parts of the variation drawn.
    lang
        ``"en"`` or ``"fr"`` (default: the global option).

    Returns
    -------
    matplotlib.figure.Figure
        The figure; save it with :func:`pyxtpanel.plot.save_figure`.

    Raises
    ------
    SpecificationError
        If ``part`` is not ``"between"``, ``"within"`` or ``"both"``.
    """
    chosen = resolve_lang(lang)
    if part not in ("between", "within", "both"):
        raise SpecificationError(
            _i18n.message("plot.scatter.bad_part", part=part), key="plot.scatter.bad_part"
        )
    (xname,), xv = term_columns(panel, x)
    (yname,), yv = term_columns(panel, y)
    parts = ("between", "within") if part == "both" else (part,)
    with chart_style():
        fig, axes = new_figure(height=3.2, shape=(1, len(parts)), squeeze=False)
        for ax, one in zip(axes.ravel(), parts, strict=True):
            data = scatter_data(panel, xv[:, 0], yv[:, 0], one)  # type: ignore[arg-type]
            colour = BLUE if one == "between" else ORANGE
            size = 14 if data.x.size < 400 else 4
            ax.scatter(
                data.x,
                data.y,
                s=size,
                color=colour,
                alpha=0.8,
                edgecolor="white",
                linewidth=0.3,
                zorder=3,
            )
            if np.isfinite(data.slope):
                grid = np.array([np.min(data.x), np.max(data.x)])
                ax.plot(
                    grid, data.intercept + data.slope * grid, color=INK, linewidth=1.1, zorder=4
                )
            ax.set_title(
                label(f"plot.scatter.{one}", chosen, name=xname), fontsize=8.5, color=MUTED
            )
            ax.set_xlabel(
                xname if one == "between" else f"{xname} ({'dev.' if chosen == 'en' else 'écart'})"
            )
            ax.set_ylabel(yname)
            ax.text(
                0.02,
                0.97,
                label("plot.scatter.fit", chosen, slope=data.slope, n=data.x.size),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7.5,
                color=MUTED,
            )
        return finish(fig, chosen)
