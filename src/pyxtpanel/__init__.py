"""pyxtpanel: panel-data management and exploration, the foundation of the pyxt estimation packages.

Importing the package loads layer 0 only; public names are loaded on first use (PEP 562).
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from pyxtpanel._version import __version__

if TYPE_CHECKING:
    from pyxtpanel._config import get_options, option_context, set_options
    from pyxtpanel._errors import (
        DataError,
        DeclarationError,
        NumericalError,
        PyxtpanelError,
        PyxtpanelWarning,
        SpecificationError,
        TimeAxisError,
    )
    from pyxtpanel.api.functions import (
        PanelLike,
        as_panel,
        attrition,
        audit,
        autocorr,
        balance_cost,
        causality,
        cointegration,
        complete_cases,
        cross_section_dependence,
        event_time,
        gap_types,
        level_shifts,
        missing_report,
        operate,
        residual_variation,
        resolve_duplicates,
        serial_correlation,
        slope_heterogeneity,
        spells,
        to_array3d,
        to_long,
        to_wide,
        transform,
        unbalancedness,
        unit_root,
        usable_equations,
        variation,
        within_outliers,
        xtcorr,
        xtdescribe,
        xtset,
        xtsum,
        xttab,
        xttrans,
    )
    from pyxtpanel.api.panel_data import PanelData
    from pyxtpanel.diagnostics import TestResult, list_tests, run_test

#: Public name -> module that defines it.
_LAZY = {
    "xtset": "pyxtpanel.api.functions",
    "resolve_duplicates": "pyxtpanel.api.functions",
    "to_long": "pyxtpanel.api.functions",
    "to_array3d": "pyxtpanel.api.functions",
    "to_wide": "pyxtpanel.api.functions",
    "level_shifts": "pyxtpanel.api.functions",
    "within_outliers": "pyxtpanel.api.functions",
    "event_time": "pyxtpanel.api.functions",
    "spells": "pyxtpanel.api.functions",
    "residual_variation": "pyxtpanel.api.functions",
    "variation": "pyxtpanel.api.functions",
    "PanelLike": "pyxtpanel.api.functions",
    "as_panel": "pyxtpanel.api.functions",
    "attrition": "pyxtpanel.api.functions",
    "audit": "pyxtpanel.api.functions",
    "autocorr": "pyxtpanel.api.functions",
    "balance_cost": "pyxtpanel.api.functions",
    "usable_equations": "pyxtpanel.api.functions",
    "causality": "pyxtpanel.api.functions",
    "cointegration": "pyxtpanel.api.functions",
    "complete_cases": "pyxtpanel.api.functions",
    "gap_types": "pyxtpanel.api.functions",
    "unbalancedness": "pyxtpanel.api.functions",
    "missing_report": "pyxtpanel.api.functions",
    "operate": "pyxtpanel.api.functions",
    "transform": "pyxtpanel.api.functions",
    "xtcorr": "pyxtpanel.api.functions",
    "xtdescribe": "pyxtpanel.api.functions",
    "xtsum": "pyxtpanel.api.functions",
    "unit_root": "pyxtpanel.api.functions",
    "cross_section_dependence": "pyxtpanel.api.functions",
    "serial_correlation": "pyxtpanel.api.functions",
    "slope_heterogeneity": "pyxtpanel.api.functions",
    "xttab": "pyxtpanel.api.functions",
    "xttrans": "pyxtpanel.api.functions",
    "duplicates": "pyxtpanel.api.functions",
    "read_stata": "pyxtpanel.api.functions",
    "read_table": "pyxtpanel.api.functions",
    "PanelData": "pyxtpanel.api.panel_data",
    "TestResult": "pyxtpanel.diagnostics",
    "list_tests": "pyxtpanel.diagnostics",
    "run_test": "pyxtpanel.diagnostics",
    "get_options": "pyxtpanel._config",
    "set_options": "pyxtpanel._config",
    "option_context": "pyxtpanel._config",
    "PyxtpanelError": "pyxtpanel._errors",
    "PyxtpanelWarning": "pyxtpanel._errors",
    "DataError": "pyxtpanel._errors",
    "DeclarationError": "pyxtpanel._errors",
    "TimeAxisError": "pyxtpanel._errors",
    "SpecificationError": "pyxtpanel._errors",
    "NumericalError": "pyxtpanel._errors",
}

__all__ = [
    "DataError",
    "DeclarationError",
    "NumericalError",
    "PanelData",
    "PanelLike",
    "PyxtpanelError",
    "PyxtpanelWarning",
    "SpecificationError",
    "TestResult",
    "TimeAxisError",
    "__version__",
    "as_panel",
    "attrition",
    "audit",
    "autocorr",
    "balance_cost",
    "causality",
    "cointegration",
    "complete_cases",
    "duplicates",
    "event_time",
    "gap_types",
    "get_options",
    "level_shifts",
    "list_tests",
    "missing_report",
    "operate",
    "option_context",
    "read_stata",
    "read_table",
    "residual_variation",
    "resolve_duplicates",
    "run_test",
    "set_options",
    "spells",
    "to_array3d",
    "to_long",
    "to_wide",
    "transform",
    "unbalancedness",
    "usable_equations",
    "variation",
    "within_outliers",
    "xtcorr",
    "xtdescribe",
    "xtset",
    "xtsum",
    "unit_root",
    "cross_section_dependence",
    "serial_correlation",
    "slope_heterogeneity",
    "xttab",
    "xttrans",
]


def __getattr__(name: str) -> Any:
    if name in _LAZY:
        value = getattr(import_module(_LAZY[name]), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'pyxtpanel' has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted([*globals(), *_LAZY])
