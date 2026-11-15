"""Metadata of the operators and transformations, and their registry (layer 3).

Every operator or transformation is a function ``f(index, values, ...)`` that takes the sorted
index of a panel and a ``float64`` array in sorted order (``(n_obs,)`` or ``(n_obs, k)``) and
returns an array of the same shape.  It registers itself with :func:`operator`; adding one is one
file in ``ops/``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal, TypeVar

from pyxtpanel._registry import Registry

F = TypeVar("F", bound=Callable[..., object])


@dataclass(frozen=True)
class OpSpec:
    """What an operator or transformation is and how it treats missing values.

    Attributes
    ----------
    kind
        ``"shift"`` (value of another period), ``"difference"`` (value minus another period),
        ``"transform"`` (uses the other observations of the unit or of the period).
    missing
        One sentence: when the result is missing.
    stata
        The Stata operator or command that is reproduced (``None`` for an addition of the
        package).
    summary
        One-line description (first line of the documentation).
    """

    kind: Literal["shift", "difference", "transform"]
    missing: str
    stata: str | None
    summary: str


#: Registry of the operators and transformations (``@operator("lag", OpSpec(...))``).
OPERATORS: Registry[OpSpec] = Registry("operator", package="pyxtpanel.ops")


def operator(name: str, spec: OpSpec) -> Callable[[F], F]:
    """Register an operator or a transformation under ``name``."""
    return OPERATORS.register(name, spec)
