"""Contract of the exploration tools (block K): function, method, metadata, documentation."""

from __future__ import annotations

import importlib
import inspect

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel.explore._spec import TOOLS

#: Arguments for each registered tool; a new tool must add its line here.
SAMPLE_CALLS: dict[str, tuple[tuple[object, ...], dict[str, object]]] = {
    "audit": ((), {"dependent": "y"}),
    "to_array3d": ((["x"],), {}),
    "to_wide": ((["x"],), {}),
    "level_shifts": (("x",), {"window": 1}),
    "within_outliers": (("x",), {"min_obs": 3}),
    "event_time": ((), {"indicator": "d"}),
    "spells": (("d",), {}),
    "residual_variation": ((["x"],), {}),
    "variation": (("x",), {}),
    "gap_types": ((), {}),
    "unbalancedness": ((), {}),
    "balance_cost": ((), {}),
    "usable_equations": (("y",), {"lags": 1}),
}


def probe() -> xt.PanelData:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 1, 2, 2, 2, 3, 3],
            "t": [1, 2, 3, 4, 1, 2, 4, 2, 3],
            "y": [1.0, 2, 4, 3, 5, 5, 7, 2, 2],
            "x": [1.0, 2, 4, 3, 5, 5, 7, 2, 2],
            "d": [0.0, 0, 1, 1, 0, 0, 1, 0, 0],
        }
    )
    return xt.xtset(frame, "id", "t")


def same(left: object, right: object) -> None:
    assert type(left) is type(right)
    if isinstance(left, pd.DataFrame):
        pd.testing.assert_frame_equal(left, right)  # type: ignore[arg-type]
        return
    if hasattr(left, "observed"):  # Array3D
        np.testing.assert_array_equal(left.values, right.values)  # type: ignore[attr-defined]
        np.testing.assert_array_equal(left.observed, right.observed)  # type: ignore[attr-defined]
        return
    for name in (
        "by_unit",
        "gaps",
        "counts",
        "windows",
        "by_row",
        "table",
        "spells",
        "coverage",
        "facts",
        "variables",
    ):
        if hasattr(left, name):
            pd.testing.assert_frame_equal(getattr(left, name), getattr(right, name))
    pd.testing.assert_frame_equal(left.to_frame(), right.to_frame())  # type: ignore[attr-defined]


@pytest.mark.parametrize("name", TOOLS.names())
def test_each_tool_is_a_function_a_method_and_gives_the_same_result(name: str) -> None:
    assert name in xt.__all__
    if not TOOLS.spec(name).on_panel:  # works on a plain frame: top-level function only
        assert not hasattr(xt.PanelData, name), f"{name} works on a frame: no PanelData method"
        return
    assert name in SAMPLE_CALLS, f"add a sample call for the new tool {name!r}"
    assert hasattr(xt.PanelData, name)
    panel = probe()
    args, kwargs = SAMPLE_CALLS[name]
    same(getattr(xt, name)(panel, *args, **kwargs), getattr(panel, name)(*args, **kwargs))


@pytest.mark.parametrize("name", TOOLS.names())
def test_each_tool_states_its_equivalents_and_how_it_is_validated(name: str) -> None:
    spec = TOOLS.spec(name)
    assert spec.summary and spec.validation
    func = TOOLS.get(name)
    assert inspect.getdoc(func)
    module = inspect.getsource(importlib.import_module(func.__module__))  # type: ignore[attr-defined]
    assert "validated" in module.lower(), f"{name}: the module must say how it is validated"
    if spec.stata:
        assert "Stata" in module
