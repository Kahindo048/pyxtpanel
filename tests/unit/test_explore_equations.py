"""Usable equations of a dynamic-panel estimator (K.2.4): by hand, loop oracle, invariances."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt
from pyxtpanel._errors import SpecificationError

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")

# Calendar 1..5.  Unit 1 complete; 2 skips t=3; 3 stops at t=3; 4 enters at t=2, y missing at t=4.
ROWS = [
    (1, 1, 1.0),
    (1, 2, 2.0),
    (1, 3, 4.0),
    (1, 4, 3.0),
    (1, 5, 6.0),
    (2, 1, 2.0),
    (2, 2, 3.0),
    (2, 4, 5.0),
    (2, 5, 4.0),
    (3, 1, 1.0),
    (3, 2, 2.0),
    (3, 3, 2.0),
    (4, 2, 3.0),
    (4, 3, 4.0),
    (4, 4, math.nan),
    (4, 5, 2.0),
]


def hand() -> xt.PanelData:
    frame = pd.DataFrame(ROWS, columns=["id", "t", "y"])
    frame["x"] = 1.0
    frame.loc[(frame.id == 1) & (frame.t == 3), "x"] = np.nan
    return xt.xtset(frame, "id", "t")


def usable(result) -> dict[int, list[int]]:  # type: ignore[no-untyped-def]
    out: dict[int, list[int]] = {}
    for (unit, period), ok in result.by_row["usable"].items():
        if ok:
            out.setdefault(int(unit), []).append(int(period))
    return out


def test_first_differences_with_one_lag_by_hand() -> None:
    result = hand().usable_equations("y", lags=1)
    # D.y_t and D.L.y_t need y at t, t-1, t-2 on the calendar
    assert usable(result) == {1: [3, 4, 5], 3: [3]}
    assert (result.n_usable, result.n_units_usable, result.n_observed) == (4, 2, 15)
    assert result.by_unit["n_usable"].to_dict() == {1: 3, 2: 0, 3: 1, 4: 0}
    assert result.by_period["usable"].tolist() == [0, 0, 2, 1, 1]
    assert result.n_instrumented == 4 and result.instrument_lag == 2


def test_forward_orthogonal_deviations_keep_rows_across_a_gap() -> None:
    result = hand().usable_equations("y", lags=1, transform="fod")
    # unit 2 is usable at t=2 although t=3 is missing: its later observations are used
    assert usable(result) == {1: [2, 3, 4], 2: [2], 3: [2], 4: [3]}
    assert (result.n_usable, result.n_units_usable) == (6, 4)
    assert result.n_instrumented == 6 and result.instrument_lag == 1


def test_a_static_model_needs_only_the_current_and_previous_period() -> None:
    result = hand().usable_equations("y", lags=0)
    assert result.n_usable == 9
    assert result.n_instrumented == 5  # a level 2 periods earlier must exist
    assert hand().usable_equations("y", lags=0, instrument_lag=1).n_instrumented == 9


def test_other_regressors_must_exist_after_the_transformation_too() -> None:
    # x is missing for unit 1 at t=3: D.x is missing at t=3 and t=4
    result = hand().usable_equations("y", lags=0, exog=["x"])
    assert result.n_usable == 7 and usable(result)[1] == [2, 5]
    assert result.exog == ("x",)


def test_the_comparison_table_covers_both_transformations() -> None:
    compare = hand().usable_equations("y", lags=1).compare
    assert compare.loc[("fd", 0), "equations"] == 9 and compare.loc[("fd", 1), "equations"] == 4
    assert compare.loc[("fod", 1), "equations"] == 6 and compare.loc[("fod", 1), "units"] == 4


def test_the_summary_and_table_forms() -> None:
    result = hand().usable_equations("y")
    assert "4 of 15" in result.summary("en") and "4 sur 15" in result.summary("fr")
    assert list(result.to_frame().columns) == ["observed", "usable", "instrumented"]
    assert "tabular" in result.to_latex()


def test_wrong_arguments_are_refused_with_the_stata_code() -> None:
    with pytest.raises(SpecificationError) as bad:
        hand().usable_equations("y", transform="gmm")  # type: ignore[arg-type]
    assert bad.value.stata_code == 198
    with pytest.raises(SpecificationError):
        hand().usable_equations("y", lags=-1)


# ----------------------------------------------------------------- independent loop oracle


def oracle(
    data: dict[int, dict[int, float]], lags: int, transform: str, m: int
) -> set[tuple[int, int]]:
    """Usable (unit, period) pairs by the definition, with dictionaries and loops only."""

    def valid(unit: int, t: int) -> bool:
        return t in data[unit] and not math.isnan(data[unit][t])

    found = set()
    for unit, rows in data.items():
        for t in rows:
            if transform == "fd":
                ok = all(valid(unit, t - k) for k in range(lags + 2))
            else:
                ok = True
                for lag_ in range(lags + 1):
                    # the series y_{t-lag} exists at row s when s is a row and y_{s-lag} is valid
                    here = valid(unit, t - lag_)
                    later = any(valid(unit, s - lag_) for s in rows if s > t)
                    ok = ok and here and later
            if ok and any(valid(unit, s) for s in rows if s <= t - m):
                found.add((unit, t))
    return found


cells = st.one_of(st.floats(-5, 5), st.just(math.nan))
panels = st.dictionaries(
    st.integers(0, 6),
    st.dictionaries(st.integers(0, 11), cells, min_size=1, max_size=12),
    min_size=1,
    max_size=7,
)


def build(data: dict[int, dict[int, float]], name=lambda u: u, shift: int = 0) -> xt.PanelData:  # type: ignore[no-untyped-def]
    rows = [(name(u), t + shift, y) for u, rs in data.items() for t, y in rs.items()]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "y"]), "id", "t")


@quiet
@settings(max_examples=150, deadline=None)
@given(data=panels, lags=st.integers(0, 3), transform=st.sampled_from(["fd", "fod"]))
def test_the_vectorised_count_equals_the_loop_implementation(data, lags, transform) -> None:  # type: ignore[no-untyped-def]
    result = build(data).usable_equations("y", lags=lags, transform=transform)
    m = 2 if transform == "fd" else 1
    expected = oracle(data, lags, transform, m)
    by_row = result.by_row
    got_usable = {(int(u), int(t)) for (u, t), ok in by_row["usable"].items() if ok}
    got_instrumented = {(int(u), int(t)) for (u, t), ok in by_row["instrumented"].items() if ok}
    assert got_instrumented == expected
    transformed = oracle(data, lags, transform, -(10**6))  # without the instrument condition
    assert got_usable == transformed
    assert result.n_usable == len(transformed) and result.n_instrumented == len(expected)


@quiet
@settings(max_examples=50, deadline=None)
@given(data=panels, shift=st.integers(-300, 300), lags=st.integers(0, 2))
def test_the_count_does_not_depend_on_unit_names_or_the_origin_of_time(data, shift, lags) -> None:  # type: ignore[no-untyped-def]
    for transform in ("fd", "fod"):
        a = build(data).usable_equations("y", lags=lags, transform=transform)  # type: ignore[arg-type]
        b = build(data, name=lambda u: f"u{u:02d}", shift=shift).usable_equations(
            "y",
            lags=lags,
            transform=transform,  # type: ignore[arg-type]
        )
        assert (a.n_usable, a.n_instrumented, a.n_units_usable) == (
            b.n_usable,
            b.n_instrumented,
            b.n_units_usable,
        )
        assert a.by_period["usable"].tolist() == b.by_period["usable"].tolist()
