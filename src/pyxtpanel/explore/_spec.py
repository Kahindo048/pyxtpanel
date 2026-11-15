"""Metadata and registry of the exploration tools (layer 4).

An exploration tool is a function ``f(panel, ...) -> Result`` registered with :func:`tool`.  Its
spec says what the nearest equivalent elsewhere is (``None`` when there is none) and how the
tool is validated, so that every tool states both in one machine-checked place (K.7.5).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TypeVar

from pyxtpanel._registry import Registry

F = TypeVar("F", bound=Callable[..., object])


@dataclass(frozen=True)
class ToolSpec:
    """What an exploration tool is and how it is validated.

    Attributes
    ----------
    summary
        One-line description.
    validation
        How the results are checked (hand-computed panels, independent implementation,
        comparison with a published formula or with another software).
    stata
        The Stata command that does the same, ``None`` when Stata has none.
    similar
        Nearest function in R or Python (``"plm::punbalancedness"``), ``None`` when none.
    on_panel
        ``True`` for a tool that takes a declared panel and has a ``PanelData`` method;
        ``False`` for one that works on a plain DataFrame (before or after a declaration).
    """

    summary: str
    validation: str
    stata: str | None = None
    similar: str | None = None
    on_panel: bool = True


#: Registry of the tools (``@tool("gap_types", ToolSpec(...))``).
TOOLS: Registry[ToolSpec] = Registry("tool", package="pyxtpanel.explore")


def tool(name: str, spec: ToolSpec) -> Callable[[F], F]:
    """Register an exploration tool under ``name``."""
    return TOOLS.register(name, spec)
