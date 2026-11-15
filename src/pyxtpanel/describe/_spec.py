"""Metadata and registry of the descriptive statistics (layer 4).

A description is a function ``f(panel, ...) -> Result`` registered with :func:`description`;
adding one is one file in ``describe/``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TypeVar

from pyxtpanel._registry import Registry

F = TypeVar("F", bound=Callable[..., object])


@dataclass(frozen=True)
class DescribeSpec:
    """What a description is and where it comes from.

    Attributes
    ----------
    stata
        The Stata command that is reproduced (``None`` for an addition of the package).
    summary
        One-line description.
    """

    stata: str | None
    summary: str


#: Registry of the descriptions (``@description("xtsum", DescribeSpec(...))``).
DESCRIPTIONS: Registry[DescribeSpec] = Registry("description", package="pyxtpanel.describe")


def description(name: str, spec: DescribeSpec) -> Callable[[F], F]:
    """Register a description under ``name``."""
    return DESCRIPTIONS.register(name, spec)
