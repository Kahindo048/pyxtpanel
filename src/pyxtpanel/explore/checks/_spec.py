"""Alerts and registry of the audit checks (layer 4).

A check is a function ``f(context) -> list[Alert]`` registered with :func:`check`.  It reads the
shared :class:`~pyxtpanel.explore.checks._context.AuditContext` (so that no result is computed
twice) and returns the alerts it raises.  Adding a check is one file in ``explore/checks/``: the
audit finds it by itself.

An alert has a stable ``code``; its text and its advice are the catalogue messages
``explore.audit.<code>`` and ``explore.audit.<code>.advice`` in French and English.  The spec of
a check lists the codes it can raise, so that a test verifies that every message exists in both
languages.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, TypeVar

from pyxtpanel import _i18n
from pyxtpanel._registry import Registry
from pyxtpanel._types import Lang

if TYPE_CHECKING:
    from pyxtpanel.explore.checks._context import AuditContext

Level = Literal["critical", "warning", "info"]
#: Levels from the most to the least serious.
LEVELS: tuple[Level, ...] = ("critical", "warning", "info")

F = TypeVar("F", bound=Callable[..., object])


@dataclass(frozen=True)
class Alert:
    """One finding of the audit.

    Attributes
    ----------
    check
        Name of the check that raised it.
    code
        Stable identifier (``few_units``, ``internal_gaps``...).
    level
        ``"critical"`` (the analysis is likely wrong or impossible), ``"warning"`` (to look at
        before estimating) or ``"info"``.
    fields
        The figures behind the message, also useful to a program.
    """

    check: str
    code: str
    level: Level
    fields: Mapping[str, object] = field(default_factory=dict)

    def message(self, lang: Lang | None = None) -> str:
        """Text of the finding in ``lang`` (default: the global ``lang`` option)."""
        return _i18n.message(f"explore.audit.{self.code}", lang, **self.fields)

    def advice(self, lang: Lang | None = None) -> str:
        """Return what to do about it, in ``lang``."""
        return _i18n.message(f"explore.audit.{self.code}.advice", lang, **self.fields)


@dataclass(frozen=True)
class CheckSpec:
    """What a check looks at.

    Attributes
    ----------
    summary
        One-line description.
    codes
        Every alert code the check can raise (each has a message and an advice in the
        catalogue).
    order
        Position of the check in the report (lower first).
    """

    summary: str
    codes: tuple[str, ...]
    order: int = 50


#: Registry of the checks (``@check("structure", CheckSpec(...))``).
CHECKS: Registry[CheckSpec] = Registry("check", package="pyxtpanel.explore.checks")


def check(name: str, spec: CheckSpec) -> Callable[[F], F]:
    """Register an audit check under ``name``."""
    return CHECKS.register(name, spec)


def short_list(names: list[str], limit: int = 8) -> str:
    """Comma-separated names, cut after ``limit`` with the number left out."""
    if len(names) <= limit:
        return ", ".join(names)
    return ", ".join(names[:limit]) + f", ... (+{len(names) - limit})"


def run(name: str, context: AuditContext) -> list[Alert]:
    """Run the check ``name`` on ``context``."""
    alerts: list[Alert] = CHECKS.get(name)(context)  # type: ignore[assignment]
    return alerts
