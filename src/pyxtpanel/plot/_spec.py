"""Registry of the figure kinds and the dispatcher ``plot(panel, kind=...)`` (F.1.4, F.2.6)."""

from __future__ import annotations

import inspect
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypeVar

from pyxtpanel import _i18n, _signature
from pyxtpanel._errors import SpecificationError
from pyxtpanel._registry import Registry
from pyxtpanel._types import Lang
from pyxtpanel.plot import _plugins

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from pyxtpanel.panel.base import PanelBase

F = TypeVar("F", bound=Callable[..., object])

_i18n.register(
    "plot",
    {
        "plot.bad_call": {
            "en": "figure {kind!r}: {reason}. Arguments of the figure: {accepted}.",
            "fr": "figure {kind!r} : {reason}. Arguments de la figure : {accepted}.",
        },
    },
)


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
    """Register a function ``f(panel, ..., lang=None) -> Figure`` as the figure kind ``name``.

    Parameters
    ----------
    name
        Name of the figure kind (``"lines"``...), the ``kind`` of :func:`pyxtpanel.plot.plot`.
    spec
        Its :class:`PlotSpec`: one-line summary and the Stata command it matches.

    Returns
    -------
    callable
        The decorator that registers the function.
    """
    return KINDS.register(name, spec)


def _resolve(kind: str, lang: Lang | None) -> Callable[..., Any]:
    """Return the function drawing ``kind``: the package's own, else one a plugin supplies."""
    if kind in KINDS:
        return KINDS.get(kind)
    draw = _plugins.find_kind(kind, lang)
    if draw is not None:
        return draw
    available = ", ".join(sorted({*KINDS.names(), *_plugins.plugin_names()}))
    key = "plot.unknown_kind" if _plugins.has_plugins() else "plot.unknown_kind_install"
    raise SpecificationError(
        _i18n.message(key, lang, kind=kind, available=available),
        key=key,
    )


def plot(panel: PanelBase, kind: str, **options: Any) -> Figure:
    """Draw the figure ``kind`` of a panel (see ``KINDS.names()`` for the kinds).

    Parameters
    ----------
    panel
        The declared panel.
    kind
        ``"lines"``, ``"scatter"``, ``"distribution"``, ``"participation"`` or
        ``"transition"``, or a kind supplied by an installed package (the entry-point group
        ``pyxtpanel.plot_kinds``, such as ``pyxtviz``); the package's own kinds come first.
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
        If the kind is unknown (the message lists the available ones and, when no package
        supplies kinds, names ``pyxtviz``), or if the call lacks an argument of the kind or gives
        one it does not accept (the message lists the arguments of the kind).
    PyxtpanelError
        If the package that supplies the kind cannot be loaded.
    """
    from pyxtpanel.plot.style import require_matplotlib

    require_matplotlib()
    lang = options.get("lang")  # the language asked for the figure also words its errors
    language: Lang | None = "en" if lang == "en" else "fr" if lang == "fr" else None
    draw = _resolve(kind, language)
    signature = inspect.signature(draw)
    try:
        signature.bind(panel, **options)
    except TypeError as error:
        raise SpecificationError(
            _i18n.message(
                "plot.bad_call",
                language,
                kind=kind,
                reason=_signature.reason(signature, error),
                accepted=_signature.accepted(signature),
            ),
            key="plot.bad_call",
        ) from None
    figure: Figure = draw(panel, **options)
    assert not isinstance(figure, type)
    return figure
