"""Invariance properties of every exploration tool (K.7.4).

Three transformations that must not change what a tool reports about a panel:

- shuffling the rows of the raw frame (the declaration sorts them);
- renaming the units (an order-preserving map, so the unit order is the same);
- shifting the time origin (the calendar step is unchanged).

One test covers all the registered tools through ``SAMPLE_CALLS`` of the contract; a tool with no
sample call fails ``test_every_panel_tool_is_covered``.  Each result is reduced to its numbers
(``canon``): the labels that carry a unit name or a period are dropped, since they move with the
transformation by construction.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt
from pyxtpanel._errors import PyxtpanelError
from pyxtpanel.explore._spec import TOOLS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "contract"))
from test_explore_tools import SAMPLE_CALLS  # noqa: E402  (shared sample calls)

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")

#: Label columns of ``to_frame()`` that name a unit or a period: they move with the transformation.
LABELS = {"unit", "time", "period", "start", "end"}

#: One unit: periods with ``(y, x, d)``; ``d`` is a 0/1 state.
cell = st.tuples(
    st.integers(-6, 6).map(float), st.integers(-6, 6).map(float), st.integers(0, 1).map(float)
)
panels = st.dictionaries(
    st.integers(0, 5),
    st.dictionaries(st.integers(0, 11), cell, min_size=1, max_size=9),
    min_size=1,
    max_size=6,
)


def build(
    units: dict[int, dict[int, tuple[float, float, float]]],
    *,
    shift: int = 0,
    rename: bool = False,
    seed: int | None = None,
) -> xt.PanelData:
    """Panel declared from the raw frame, optionally renamed, shifted and shuffled."""
    rows = [
        (f"unit_{u:02d}" if rename else u, t + shift, *vals)
        for u, periods in units.items()
        for t, vals in periods.items()
    ]
    frame = pd.DataFrame(rows, columns=["id", "t", "y", "x", "d"])
    if seed is not None:
        frame = frame.iloc[np.random.default_rng(seed).permutation(len(frame))]
    return xt.xtset(frame.reset_index(drop=True), "id", "t")


def canon(result: object) -> object:
    """The numbers of a result, without unit or period labels."""
    if hasattr(result, "observed"):  # Array3D
        return (np.asarray(result.values), np.asarray(result.observed))  # type: ignore[attr-defined]
    if isinstance(result, pd.DataFrame):  # to_wide: the columns carry the period
        return np.asarray(result.to_numpy(dtype=float))
    frame = result.to_frame()  # type: ignore[attr-defined]
    if frame.index.name in LABELS or frame.index.names[0] in LABELS:
        frame = frame.reset_index(drop=True)
    if "check" in frame.columns:  # AuditReport: the messages may quote labels
        return frame[["level", "check", "code"]].reset_index(drop=True)
    return frame.drop(columns=[c for c in frame.columns if c in LABELS]).reset_index(drop=True)


def attempt(panel: xt.PanelData, name: str) -> tuple[str, object]:
    args, kwargs = SAMPLE_CALLS[name]
    try:
        return "ok", canon(getattr(panel, name)(*args, **kwargs))
    except PyxtpanelError as error:
        return "error", error.code if hasattr(error, "code") else type(error).__name__


def same(left: tuple[str, object], right: tuple[str, object]) -> None:
    assert left[0] == right[0]
    a, b = left[1], right[1]
    if isinstance(a, pd.DataFrame):
        pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=1e-10)  # type: ignore[arg-type]
    elif isinstance(a, tuple):
        for x, y in zip(a, b, strict=True):  # type: ignore[arg-type]
            np.testing.assert_allclose(x, y, rtol=1e-10, equal_nan=True)
    elif isinstance(a, np.ndarray):
        np.testing.assert_allclose(a, b, rtol=1e-10, equal_nan=True)  # type: ignore[arg-type]
    else:
        assert a == b


PANEL_TOOLS = [name for name in TOOLS.names() if TOOLS.spec(name).on_panel]


def test_every_panel_tool_is_covered() -> None:
    assert set(PANEL_TOOLS) <= set(SAMPLE_CALLS)


@quiet
@pytest.mark.parametrize("name", PANEL_TOOLS)
@settings(max_examples=40, deadline=None)
@given(units=panels, seed=st.integers(0, 10_000))
def test_shuffling_the_rows_changes_nothing(name: str, units, seed) -> None:  # type: ignore[no-untyped-def]
    same(attempt(build(units), name), attempt(build(units, seed=seed), name))


@quiet
@pytest.mark.parametrize("name", PANEL_TOOLS)
@settings(max_examples=40, deadline=None)
@given(units=panels)
def test_renaming_the_units_changes_nothing(name: str, units) -> None:  # type: ignore[no-untyped-def]
    same(attempt(build(units), name), attempt(build(units, rename=True), name))


@quiet
@pytest.mark.parametrize("name", PANEL_TOOLS)
@settings(max_examples=40, deadline=None)
@given(units=panels, shift=st.integers(-500, 500))
def test_shifting_the_time_changes_nothing(name: str, units, shift) -> None:  # type: ignore[no-untyped-def]
    same(attempt(build(units), name), attempt(build(units, shift=shift), name))


# ------------------------------------------------------------- tools on a plain frame

duplicated = st.lists(
    st.tuples(st.integers(0, 3), st.integers(0, 3), st.integers(-5, 5).map(float)),
    min_size=1,
    max_size=14,
)


@quiet
@settings(max_examples=60, deadline=None)
@given(
    rows=duplicated,
    how=st.sampled_from(["mean", "sum", "min", "max", "median"]),
    seed=st.integers(0, 10_000),
    shift=st.integers(-100, 100),
)
def test_resolving_duplicates_by_a_statistic_ignores_order_names_and_origin(  # type: ignore[no-untyped-def]
    rows, how, seed, shift
) -> None:
    frame = pd.DataFrame(rows, columns=["id", "t", "x"])
    base = xt.resolve_duplicates(frame, "id", "t", how=how)
    moved = frame.assign(id=frame["id"].map(lambda u: f"unit_{u:02d}"), t=frame["t"] + shift)
    moved = moved.iloc[np.random.default_rng(seed).permutation(len(moved))]
    other = xt.resolve_duplicates(moved, "id", "t", how=how)
    key = ["id", "t"]
    a = base.frame.sort_values(key).reset_index(drop=True)
    b = other.frame.sort_values(key).reset_index(drop=True)
    np.testing.assert_allclose(a["x"], b["x"], rtol=1e-12, equal_nan=True)
    assert (a["t"] + shift).tolist() == b["t"].tolist()
    assert base.n_rows_in == other.n_rows_in and base.n_rows_out == other.n_rows_out


@quiet
@settings(max_examples=60, deadline=None)
@given(units=panels, seed=st.integers(0, 10_000), shift=st.integers(-100, 100))
def test_wide_and_long_ignore_the_row_order_and_the_origin_of_time(units, seed, shift) -> None:  # type: ignore[no-untyped-def]
    base = build(units).to_wide(["x"]).reset_index()
    other = build(units, seed=seed, shift=shift).to_wide(["x"]).reset_index()
    np.testing.assert_array_equal(
        base.drop(columns="id").to_numpy(float), other.drop(columns="id").to_numpy(float)
    )
    long_base = xt.to_long(base, ["x"], id="id", time="t", drop_empty=True)
    long_other = xt.to_long(other, ["x"], id="id", time="t", drop_empty=True)
    assert (long_base["t"] + shift).tolist() == long_other["t"].tolist()
    np.testing.assert_array_equal(long_base["x"].to_numpy(), long_other["x"].to_numpy())
