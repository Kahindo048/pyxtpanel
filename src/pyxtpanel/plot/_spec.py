"""Registry of the figure kinds and the dispatcher ``plot(panel, kind=...)`` (F.1.4, F.2.6)."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypeVar

from pyxtpanel._registry import Registry

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from pyxtpanel.panel.base import PanelBase

F = TypeVar("F", bound=Callable[..., object])


@dataclass(frozen=True)
class PlotSpec:
    """What a kind of figure shows.

    Attributes
    ----------
    summary
        One-line description.
    stata
        The Stata command it corresponds to (``None`` for an addition of the package).
    """

    summary: str
    stata: str | None = None


#: Registry of the figure kinds (``@plot_kind("lines", PlotSpec(...))``).
KINDS: Registry[PlotSpec] = Registry("figure kind", package="pyxtpanel.plot.kinds")


def plot_kind(name: str, spec: PlotSpec) -> Callable[[F], F]:
    """Register a function ``f(panel, ..., lang=None) -> Figure`` as the figure kind ``name``."""
    return KINDS.register(name, spec)


def plot(panel: PanelBase, kind: str, **options: Any) -> Figure:
    """Draw the figure ``kind`` of a panel (see ``KINDS.names()`` for the kinds).

    Parameters
    ----------
    panel
        The declared panel.
    kind
        ``"lines"``, ``"scatter"``, ``"distribution"``, ``"participation"`` or
        ``"transition"``.
    **options
        Options of the kind (``variable``, ``x`` and ``y``, ``lang``...).

    Returns
    -------
    matplotlib.figure.Figure
        The figure; save it with :func:`pyxtpanel.plot.style.save_figure`.

    Raises
    ------
    ImportError
        If matplotlib is not installed (``pip install pyxtpanel[plot]``).
    SpecificationError
        If the kind is unknown (the message lists the available ones).
    """
    from pyxtpanel.plot.style import require_matplotlib

    require_matplotlib()
    figure = KINDS.get(kind)(panel, **options)
    assert not isinstance(figure, type)
    return figure  # type: ignore[return-value]
