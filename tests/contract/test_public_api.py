"""Contract tests of the public surface: snapshot, registries, Stata provenance, error codes."""

import ast
import importlib
import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from tools import snapshot_api, stata_map

import pyxtpanel as xt
from pyxtpanel._errors import STATA_CODES
from pyxtpanel.describe._spec import DESCRIPTIONS
from pyxtpanel.ops._spec import OPERATORS

SRC = Path(xt.__file__).parent


def test_the_public_api_matches_its_snapshot() -> None:
    # a deliberate change: python -m tools.snapshot_api --write (and note it in CHANGELOG.md)
    assert snapshot_api.main([]) == 0


def test_the_stata_map_is_up_to_date() -> None:
    assert stata_map.main([]) == 0


def test_every_public_name_resolves_and_is_documented() -> None:
    for name in xt.__all__:
        value = getattr(xt, name)
        if inspect.isfunction(value) or inspect.isclass(value):
            assert inspect.getdoc(value), f"{name} has no docstring"


# ---------------------------------------------------------------- registries


def probe() -> xt.PanelData:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 1, 2, 2, 2, 3, 3],
            "t": [1, 2, 3, 4, 1, 2, 4, 2, 3],
            "x": [1.0, 2, 3, 4, 5, 5, 7, 2, 2],
            "z": [10.0, 11, 12, 13, np.nan, 15, 17, 20, 21],
            "d": [0.0, 0, 1, 1, 1, 1, 2, 0, 0],
        }
    )
    return xt.xtset(frame, "id", "t")


#: Arguments for each registered description; a new description must add its line here.
SAMPLE_CALLS = {
    "xtdescribe": ((), {}),
    "xtsum": ((["x", "z"],), {}),
    "xttab": (("d",), {}),
    "xttrans": (("d",), {}),
    "xtcorr": ((["x", "z"],), {}),
    "autocorr": (("x",), {"lags": [1, 2]}),
    "attrition": ((), {}),
    "missing_report": ((), {}),
    "complete_cases": ((["z"],), {}),
}


def _same(left: object, right: object) -> bool:
    if hasattr(left, "to_frame"):
        left, right = left.to_frame(), right.to_frame()  # type: ignore[attr-defined]
    if isinstance(left, pd.DataFrame):
        pd.testing.assert_frame_equal(left, right)  # type: ignore[arg-type]
        return True
    return bool(left == right)


@pytest.mark.parametrize("name", DESCRIPTIONS.names())
def test_each_description_is_a_function_a_method_and_gives_the_same_result(name: str) -> None:
    assert name in SAMPLE_CALLS, f"add a sample call for the new description {name!r}"
    assert name in xt.__all__ and hasattr(xt.PanelData, name)
    panel = probe()
    args, kwargs = SAMPLE_CALLS[name]
    if name == "complete_cases":  # the result holds a panel and a mask, not only a table
        a, b = getattr(xt, name)(panel, *args, **kwargs), getattr(panel, name)(*args, **kwargs)
        assert a.n_rows_kept == b.n_rows_kept and a.mask.equals(b.mask)
        return
    assert _same(getattr(xt, name)(panel, *args, **kwargs), getattr(panel, name)(*args, **kwargs))


@pytest.mark.parametrize("name", OPERATORS.names())
def test_each_operator_has_metadata_and_a_stata_or_package_origin(name: str) -> None:
    spec = OPERATORS.spec(name)
    assert spec.summary and spec.missing and spec.kind in {"shift", "difference", "transform"}
    assert inspect.getdoc(OPERATORS.get(name))


def _defining_source(func: object) -> str:
    return inspect.getsource(importlib.import_module(func.__module__))  # type: ignore[attr-defined]


@pytest.mark.parametrize("name", DESCRIPTIONS.names())
def test_a_description_that_reproduces_stata_says_so_in_its_module(name: str) -> None:
    spec = DESCRIPTIONS.spec(name)
    if spec.stata:
        assert "Stata" in _defining_source(DESCRIPTIONS.get(name))


@pytest.mark.parametrize("name", OPERATORS.names())
def test_an_operator_that_reproduces_stata_says_so_in_its_module(name: str) -> None:
    spec = OPERATORS.spec(name)
    if spec.stata:
        assert "Stata" in _defining_source(OPERATORS.get(name))


def test_the_provenance_comments_are_collected() -> None:
    found = stata_map.quirks()
    assert any(marker == "STATA-QUIRK" for _, marker, _ in found)
    assert any("xtsum" in text for _, marker, text in found if marker == "STATA-QUIRK")


# ---------------------------------------------------------------- Stata return codes


def test_every_stata_code_in_the_source_is_in_the_table() -> None:
    used: dict[int, str] = {}
    for path in SRC.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            holder = None
            if (
                isinstance(node, ast.keyword)
                and node.arg == "stata_code"
                or isinstance(node, ast.Assign)
                and any(getattr(target, "id", "") == "stata_code" for target in node.targets)
            ):
                holder = node.value
            if isinstance(holder, ast.Constant) and isinstance(holder.value, int):
                used[holder.value] = path.name
    unlisted = {code: where for code, where in used.items() if code not in STATA_CODES}
    assert not unlisted, f"codes missing from STATA_CODES: {unlisted}"
    assert {198, 451, 459, 2000} <= set(used)
