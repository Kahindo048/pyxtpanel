"""Spells and event time (K.3.2, K.3.4): by hand, loop oracle, invariances."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt
from pyxtpanel._errors import DataError, SpecificationError

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")
nan = math.nan

# ----------------------------------------------------------------- spells, by hand


def spell_panel() -> xt.PanelData:
    rows = []
    for t, d in zip(range(1, 7), [0, 0, 1, 1, 1, 0], strict=True):
        rows.append((1, t, float(d)))
    for t, d in zip([1, 2, 4, 5], [1, 1, 1, 2], strict=True):  # unit 2 skips t=3
        rows.append((2, t, float(d)))
    for t, d in zip([1, 2, 3], [5, nan, 5], strict=True):  # unit 3 has a missing value at t=2
        rows.append((3, t, float(d)))
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "d"]), "id", "t")


def test_spells_by_hand() -> None:
    result = spell_panel().spells("d")
    got = result.spells[["unit", "state", "start", "end", "length", "start_reason", "end_reason"]]
    assert got.to_dict("records") == [
        {
            "unit": 1,
            "state": 0,
            "start": 1,
            "end": 2,
            "length": 2,
            "start_reason": "first",
            "end_reason": "change",
        },
        {
            "unit": 1,
            "state": 1,
            "start": 3,
            "end": 5,
            "length": 3,
            "start_reason": "change",
            "end_reason": "change",
        },
        {
            "unit": 1,
            "state": 0,
            "start": 6,
            "end": 6,
            "length": 1,
            "start_reason": "change",
            "end_reason": "last",
        },
        {
            "unit": 2,
            "state": 1,
            "start": 1,
            "end": 2,
            "length": 2,
            "start_reason": "first",
            "end_reason": "gap",
        },
        {
            "unit": 2,
            "state": 1,
            "start": 4,
            "end": 4,
            "length": 1,
            "start_reason": "gap",
            "end_reason": "change",
        },
        {
            "unit": 2,
            "state": 2,
            "start": 5,
            "end": 5,
            "length": 1,
            "start_reason": "change",
            "end_reason": "last",
        },
        {
            "unit": 3,
            "state": 5,
            "start": 1,
            "end": 1,
            "length": 1,
            "start_reason": "first",
            "end_reason": "gap",
        },
        {
            "unit": 3,
            "state": 5,
            "start": 3,
            "end": 3,
            "length": 1,
            "start_reason": "gap",
            "end_reason": "last",
        },
    ]
    assert result.by_unit["n_changes"].to_dict() == {1: 2, 2: 1, 3: 0}
    assert result.by_unit["n_spells"].to_dict() == {1: 3, 2: 3, 3: 2}
    assert (result.n_changes, result.n_units_changing, result.n_units) == (3, 2, 3)
    # complete spells have neither end open: unit 1 state 1 (3), unit 2 state 1 (1)
    assert result.by_state.loc[1, "complete"] == 2 and result.by_state.loc[1, "mean"] == 2.0
    assert result.by_state.loc[0, "complete"] == 0 and math.isnan(result.by_state.loc[0, "mean"])
    assert "3 changes" in result.summary("en") and "3 changements" in result.summary("fr")


def test_spells_refuse_a_continuous_variable() -> None:
    frame = pd.DataFrame({"id": [1] * 4, "t": [1, 2, 3, 4], "x": [0.1, 0.2, 0.3, 0.4]})
    with pytest.raises(SpecificationError) as bad:
        xt.xtset(frame, "id", "t").spells("x", max_states=3)
    assert bad.value.stata_code == 134


def loop_spells(data: dict[int, dict[int, float]]) -> list[tuple[object, ...]]:
    out = []
    for unit, rows in data.items():
        times = sorted(rows)
        run: list[int] = []
        for t in times:
            valid = not math.isnan(rows[t])
            if valid and run and t == run[-1] + 1 and rows[t] == rows[run[0]]:
                run.append(t)
                continue
            if run:
                out.append((unit, rows[run[0]], run[0], run[-1]))
            run = [t] if valid else []
        if run:
            out.append((unit, rows[run[0]], run[0], run[-1]))
    return out


cells = st.one_of(st.integers(0, 2).map(float), st.just(nan))
panels = st.dictionaries(
    st.integers(0, 5),
    st.dictionaries(st.integers(0, 9), cells, min_size=1, max_size=10),
    min_size=1,
    max_size=6,
)


def build(data: dict[int, dict[int, float]], name=lambda u: u, shift: int = 0) -> xt.PanelData:  # type: ignore[no-untyped-def]
    rows = [(name(u), t + shift, y) for u, rs in data.items() for t, y in rs.items()]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "d"]), "id", "t")


@quiet
@settings(max_examples=150, deadline=None)
@given(data=panels)
def test_spells_equal_the_loop_implementation(data) -> None:  # type: ignore[no-untyped-def]
    if all(math.isnan(y) for rows in data.values() for y in rows.values()):
        return
    table = build(data).spells("d").spells
    got = list(zip(table["unit"], table["state"], table["start"], table["end"], strict=True))
    assert [(int(u), float(s), int(a), int(b)) for u, s, a, b in got] == [
        (u, float(s), a, b) for u, s, a, b in sorted(loop_spells(data), key=lambda r: (r[0], r[2]))
    ]


@quiet
@settings(max_examples=50, deadline=None)
@given(data=panels, shift=st.integers(-200, 200))
def test_spells_do_not_depend_on_unit_names_or_the_origin_of_time(data, shift) -> None:  # type: ignore[no-untyped-def]
    if all(math.isnan(y) for rows in data.values() for y in rows.values()):
        return
    a = build(data).spells("d")
    b = build(data, name=lambda u: f"u{u:02d}", shift=shift).spells("d")
    assert a.spells["length"].tolist() == b.spells["length"].tolist()
    assert (a.spells["start"] + shift).tolist() == b.spells["start"].tolist()
    assert a.by_unit["n_changes"].tolist() == b.by_unit["n_changes"].tolist()


# ----------------------------------------------------------------- event time, by hand


def event_panel() -> xt.PanelData:
    rows = []
    for t, d in zip(range(1, 7), [0, 0, 1, 1, 1, 1], strict=True):
        rows.append(("A", t, float(d), 3.0))
    for t, d in zip(range(2, 7), [0, 0, 0, 1, 1], strict=True):
        rows.append(("B", t, float(d), 5.0))
    for t in range(1, 7):
        rows.append(("C", t, 0.0, nan))
    for t, d in zip(range(1, 5), [0, 0, 1, 1], strict=True):
        rows.append(("D", t, float(d), 3.0))
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "d", "first"]), "id", "t")


@pytest.mark.parametrize("how", ["indicator", "event_period"])
def test_event_time_by_hand(how: str) -> None:
    panel = event_panel()
    result = panel.event_time(
        **({"indicator": "d"} if how == "indicator" else {"event_period": "first"})
    )
    unit = result.by_unit
    assert unit["treated"].to_dict() == {"A": True, "B": True, "C": False, "D": True}
    assert unit["event"].dropna().astype(int).to_dict() == {"A": 3, "B": 5, "D": 3}
    assert unit["min_rel"].dropna().to_dict() == {"A": -2, "B": -3, "D": -2}
    assert unit["max_rel"].dropna().to_dict() == {"A": 3, "B": 1, "D": 1}
    assert unit.loc[["A", "B", "D"], "pre"].tolist() == [2, 3, 2]
    assert unit.loc[["A", "B", "D"], "post"].tolist() == [4, 2, 2]
    assert result.coverage["units"].to_dict() == {-3: 1, -2: 3, -1: 3, 0: 3, 1: 3, 2: 1, 3: 1}
    assert result.coverage.loc[-3, "share"] == pytest.approx(1 / 3)
    assert (result.common_pre, result.common_post) == (2, 1)
    assert (result.n_treated, result.n_never, result.n_units) == (3, 1, 4)
    relative = result.relative
    assert relative.index.equals(panel.frame.index)
    assert relative[(panel.frame["id"] == "A").to_numpy()].tolist() == [-2, -1, 0, 1, 2, 3]
    assert relative[(panel.frame["id"] == "C").to_numpy()].isna().all()
    assert "3 of 4" in result.summary("en") and "3 unités sur 4" in result.summary("fr")


def test_event_time_needs_exactly_one_definition() -> None:
    panel = event_panel()
    for kwargs in ({}, {"indicator": "d", "event_period": "first"}):
        with pytest.raises(SpecificationError) as bad:
            panel.event_time(**kwargs)
        assert bad.value.stata_code == 198


def test_an_event_period_must_be_constant_within_the_unit() -> None:
    frame = event_panel().to_frame().reset_index(drop=True)
    frame.loc[(frame["id"] == "A") & (frame["t"] == 4), "first"] = 4.0
    with pytest.raises(DataError) as bad:
        xt.xtset(frame, "id", "t").event_time(event_period="first")
    assert bad.value.stata_code == 459


def test_event_time_on_a_monthly_calendar() -> None:
    months = pd.date_range("2020-01-01", periods=5, freq="MS")
    frame = pd.DataFrame({"id": 1, "t": months, "when": pd.Timestamp("2020-03-01")})
    result = xt.xtset(frame, "id", "t", freq="M").event_time(event_period="when")
    assert result.relative.tolist() == [-2, -1, 0, 1, 2]
    assert result.by_unit.loc[1, "event"] == pd.Timestamp("2020-03-01")


def test_with_a_calendar_step_the_relative_time_counts_steps() -> None:
    frame = pd.DataFrame({"id": [1] * 4, "t": [0, 5, 10, 15], "d": [0.0, 0, 1, 1]})
    result = xt.xtset(frame, "id", "t", delta=5).event_time(indicator="d")
    assert result.relative.tolist() == [-2, -1, 0, 1]


def loop_relative(data: dict[int, dict[int, float]]) -> dict[tuple[int, int], float]:
    out = {}
    for unit, rows in data.items():
        hits = [t for t, y in rows.items() if y != 0 and not math.isnan(y)]
        for t in rows:
            out[(unit, t)] = (t - min(hits)) if hits else nan
    return out


event_cells = st.sampled_from([0.0, 0.0, 1.0, nan])
event_panels = st.dictionaries(
    st.integers(0, 5),
    st.dictionaries(st.integers(0, 9), event_cells, min_size=1, max_size=10),
    min_size=1,
    max_size=6,
)


@quiet
@settings(max_examples=120, deadline=None)
@given(data=event_panels)
def test_event_time_equals_the_loop_implementation(data) -> None:  # type: ignore[no-untyped-def]
    rows = [(u, t, y) for u, rs in data.items() for t, y in rs.items()]
    frame = pd.DataFrame(rows, columns=["id", "t", "d"])
    result = xt.xtset(frame, "id", "t").event_time(indicator="d")
    expected = loop_relative(data)
    got = dict(
        zip(zip(frame["id"], frame["t"], strict=True), result.relative.to_numpy(), strict=True)
    )
    for key, value in expected.items():
        np.testing.assert_equal(got[key], value)
    treated = [u for u in data if any(y != 0 and not math.isnan(y) for y in data[u].values())]
    assert result.n_treated == len(treated)
