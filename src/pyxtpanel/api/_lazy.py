"""Deferred import of the statistical tests (layer 5).

The tests need ``scipy.stats``, which takes about a second to import.  The API modules reach the
tests through :func:`run_panel_test`, which imports them at the first call, so that declaring a
panel or describing it does not pay for them.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyxtpanel.diagnostics import TestResult


def run_panel_test(
    panel: object,
    test: str | None,
    *args: object,
    family: str | None = None,
    **options: object,
) -> TestResult:
    """Run a category A test (see :func:`pyxtpanel.diagnostics.run_panel_test`)."""
    from pyxtpanel.diagnostics import run_panel_test as run

    return run(panel, test, *args, family=family, **options)
