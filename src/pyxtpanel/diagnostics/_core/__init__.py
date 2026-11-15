"""TestResult, TestSpec, registry and input normalisation shared by every test.

Layer 4.
"""

from pyxtpanel.diagnostics._core.dispatch import panel_tests, run_panel_test
from pyxtpanel.diagnostics._core.equations import Equations
from pyxtpanel.diagnostics._core.inputs import as_array, panel_variables, plain_label
from pyxtpanel.diagnostics._core.result import Note, TestResult
from pyxtpanel.diagnostics._core.spec import (
    FAMILIES,
    TESTS,
    TestSpec,
    list_tests,
    register_test,
    resolve_test,
    run_test,
)

__all__ = [
    "FAMILIES",
    "Equations",
    "TESTS",
    "Note",
    "TestResult",
    "TestSpec",
    "as_array",
    "list_tests",
    "panel_tests",
    "panel_variables",
    "plain_label",
    "register_test",
    "resolve_test",
    "run_panel_test",
    "run_test",
]
