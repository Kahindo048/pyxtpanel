"""``plot(kind="distribution")``: between and within distributions of a variable (F.2.3)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._types import FloatArray, Lang
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
        "plot.dist.between": {
            "en": "Between: one point per unit (unit means)",
            "fr": "Entre : un point par unité (moyennes par unité)",
        },
        "plot.dist.within": {
            "en": "Within: x − unit mean + overall mean",
            "fr": "Intra : x − moyenne de l'unité + moyenne globale",
        },
        "plot.dist.units": {"en": "Units", "fr": "Unités"},
        "plot.dist.observations": {"en": "Observations", "fr": "Observations"},
        "plot.dist.sd": {"en": "s.d. {sd:.3g}", "fr": "écart type {sd:.3g}"},
        "plot.dist.bad_style": {
            "en": "style must be 'hist' or 'box', got {style!r}.",
            "fr": "style doit valoir 'hist' ou 'box', reçu {style!r}.",
        },
    },
)


@dataclass(frozen=True)
class DistributionData:
    """The two samples of the figure: unit means and within values."""

    between: FloatArray
    within: FloatArray


def distribution_data(panel: PanelBase, values: FloatArray) -> DistributionData:
    """Return the unit means and the within values ``x − x̄_i + x̄`` (as in ``xtsum``)."""
    keep = ~np.isnan(values)
    x, unit = values[keep], panel.index.unit[keep]
    size = np.bincount(unit, minlength=panel.n_units)
    mean = np.bincount(unit, weights=x, minlength=panel.n_units) / np.maximum(size, 1)
    return DistributionData(mean[size > 0], x - mean[unit] + x.mean())


@plot_kind(
    "distribution",
    PlotSpec(summary="Distribution of a variable between units and within units."),
)
def plot_distribution(
    panel: PanelBase,
    variable: str,
    *,
    style: Literal["hist", "box"] = "hist",
    bins: int = 30,
    lang: Lang | None = None,
) -> Figure:
    """Histograms (or boxes) of the unit means and of the within values of a variable.

    The within values are shifted by the overall mean, as in ``xtsum``, so that both panels share
    the scale of the variable.  Each panel says its standard deviation.
    """
    chosen = resolve_lang(lang)
    if style not in ("hist", "box"):
        raise SpecificationError(
            _i18n.message("plot.dist.bad_style", style=style), key="plot.dist.bad_style"
        )
    (name,), values = term_columns(panel, variable)
    data = distribution_data(panel, values[:, 0])
    panels = (
        ("between", data.between, BLUE, "plot.dist.units"),
        ("within", data.within, ORANGE, "plot.dist.observations"),
    )
    with chart_style():
        fig, axes = new_figure(height=3.0, shape=(1, 2), sharex=(style == "hist"))
        for ax, (part, sample, colour, ylabel) in zip(axes, panels, strict=True):
            if style == "hist":
                ax.hist(
                    sample, bins=bins, color=colour, alpha=0.85, edgecolor="white", linewidth=0.5
                )
                ax.set_ylabel(label(ylabel, chosen))
                ax.set_xlabel(name)
                ax.grid(axis="x", visible=False)
            else:
                options = {
                    "patch_artist": True,
                    "boxprops": {"facecolor": colour, "alpha": 0.6, "edgecolor": INK},
                    "medianprops": {"color": INK},
                    "flierprops": {"markersize": 2.5, "markeredgecolor": MUTED},
                }
                try:  # matplotlib >= 3.10 names the direction; older versions use vert=False
                    ax.boxplot(sample, orientation="horizontal", **options)
                except TypeError:
                    ax.boxplot(sample, vert=False, **options)
                ax.set_yticks([])
                ax.set_xlabel(name)
            ax.set_title(label(f"plot.dist.{part}", chosen), fontsize=8.5, color=MUTED)
            sd = float(np.std(sample, ddof=1)) if sample.size > 1 else float("nan")
            ax.text(
                0.98,
                0.96,
                label("plot.dist.sd", chosen, sd=sd),
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7.5,
                color=MUTED,
            )
        return finish(fig, chosen)
