"""Operators L F D S: parser, calendar semantics and edge cases (D.1)."""

import numpy as np
import pandas as pd
import pytest

from pyxtpanel import xtset
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.ops._spec import OPERATORS
from pyxtpanel.ops.parse import Term, parse_term, parse_terms

NAN = np.nan


def same(actual, expected):  # type: ignore[no-untyped-def]
    left, right = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    np.testing.assert_array_equal(left, right)


def make(times: list[int], x: list[float], ids: list[int] | None = None, **opts):  # type: ignore[no-untyped-def]
    ids = ids if ids is not None else [1] * len(times)
    return xtset(pd.DataFrame({"id": ids, "t": times, "x": x}), "id", "t", **opts)


# ---------------------------------------------------------------- parser


@pytest.mark.parametrize(
    ("text", "canonical"),
    [
        ("x", "x"),
        ("L.x", "L.x"),
        ("l2.x", "L2.x"),
        ("ll.x", "L2.x"),
        ("L0.x", "x"),
        ("lf.x", "x"),
        ("f2.x", "F2.x"),
        ("dl.x", "LD.x"),
        ("LD.x", "LD.x"),
        ("d2.x", "D2.x"),
        ("S12.x", "S12.x"),
        ("S.x", "S.x"),
        ("ds2.x", "DS2.x"),
    ],
)
def test_canonical_spelling(text: str, canonical: str) -> None:
    (term,) = parse_term(text)
    assert term.name == canonical


def test_ranges_expand_in_order() -> None:
    assert [t.name for t in parse_term("L(1/3).x")] == ["L.x", "L2.x", "L3.x"]
    assert [t.name for t in parse_term("L(1 3).x")] == ["L.x", "L3.x"]
    assert [t.name for t in parse_term("L(1(2)5).x")] == ["L.x", "L3.x", "L5.x"]
    assert [t.name for t in parse_terms(["L(1/2).x", "D.y"])] == ["L.x", "L2.x", "D.y"]


@pytest.mark.parametrize("bad", ["", "L.", ".x", "1x", "L(a).x", "L().x", "Q.x", "L-1.x"])
def test_invalid_terms_are_refused_with_code_198(bad: str) -> None:
    with pytest.raises(SpecificationError) as info:
        parse_term(bad)
    assert info.value.stata_code == 198 and info.value.key == "ops.term_invalid"


def test_term_properties() -> None:
    term = Term("x", shift=-1)
    assert term.has_dsf and term.netlag == -1 and term.name == "F.x"
    assert not Term("x", shift=2).has_dsf


# ---------------------------------------------------------------- calendar semantics


def test_lag_is_the_value_of_the_previous_period_not_row() -> None:
    panel = make([1, 2, 4, 5], [10.0, 20, 40, 50])
    same(panel.lag("x"), [NAN, 10.0, NAN, 40.0])
    same(panel.lag("x", 2), [NAN, NAN, 20.0, NAN])
    same(panel.lead("x"), [20.0, NAN, 50.0, NAN])


def test_diff_and_seasonal_need_both_periods() -> None:
    panel = make([1, 2, 4, 5], [10.0, 20, 40, 50])
    same(panel.diff("x"), [NAN, 10.0, NAN, 10.0])
    same(panel.diff("x", 2), [NAN, NAN, NAN, NAN])
    same(panel.seasonal("x", 2), [NAN, NAN, 20.0, NAN])
    gap = make([1, 2, 3, 5], [10.0, 20, 30, 50])
    same(gap.operate("LD.x")["LD.x"], [NAN, NAN, 10.0, NAN])


def test_second_difference_on_a_complete_unit() -> None:
    panel = make([1, 2, 3, 4], [1.0, 4, 9, 16])
    same(panel.diff("x", 2), [NAN, NAN, 2.0, 2.0])


def test_steps_dates_and_units_never_mix() -> None:
    step2 = make([1, 3, 5], [1.0, 2, 3], delta=2)
    same(step2.lag("x"), [NAN, 1.0, 2.0])
    dates = pd.DataFrame(
        {
            "id": [1, 1, 1, 2],
            "t": pd.to_datetime(["2020-01-01", "2021-01-01", "2023-01-01", "2020-01-01"]),
            "x": [1.0, 2, 3, 9],
        }
    )
    yearly = xtset(dates, "id", "t", freq="Y")
    same(yearly.lag("x"), [NAN, 1.0, NAN, NAN])
    two_units = make([1, 2, 1, 2], [1.0, 2, 5, 6], ids=[1, 1, 2, 2])
    same(two_units.lag("x"), [NAN, 1.0, NAN, 5.0])  # the last of unit 1 is not a lag


def test_results_are_aligned_on_the_original_rows_and_index() -> None:
    frame = pd.DataFrame(
        {"id": [2, 1, 1, 2], "t": [2, 2, 1, 1], "x": [4.0, 2, 1, 3]}, index=list("abcd")
    )
    panel = xtset(frame, "id", "t")
    out = panel.lag("x")
    assert out.index.tolist() == list("abcd")
    same(out, [3.0, 1.0, NAN, NAN])
    assert out.name == "L.x"
    assert frame["x"].tolist() == [4.0, 2, 1, 3]  # the user's frame is untouched


# ---------------------------------------------------------------- edge cases (D.1.4)


def test_lag_larger_than_t_and_single_observation() -> None:
    panel = make([1, 2, 3], [1.0, 2, 3])
    assert panel.lag("x", 10).isna().all() and panel.lead("x", 3).isna().all()
    single = make([5], [7.0])
    assert single.lag("x").isna().all() and single.diff("x").isna().all()


def test_negative_lag_is_a_lead_and_zero_is_the_variable() -> None:
    panel = make([1, 2, 3], [1.0, 2, 3])
    same(panel.lag("x", -1), panel.lead("x"))
    assert panel.lag("x", 0).tolist() == [1.0, 2.0, 3.0]
    with pytest.raises(SpecificationError, match="non-negative"):
        panel.diff("x", -1)


def test_missing_values_propagate() -> None:
    panel = make([1, 2, 3, 4], [1.0, NAN, 3, 4])
    same(panel.lag("x"), [NAN, 1.0, NAN, 3.0])
    same(panel.diff("x"), [NAN, NAN, NAN, 1.0])


def test_several_variables_and_amounts_give_a_frame() -> None:
    panel = xtset(
        pd.DataFrame({"id": [1, 1, 1], "t": [1, 2, 3], "x": [1.0, 2, 3], "y": [4.0, 5, 6]}),
        "id",
        "t",
    )
    out = panel.lag(["x", "y"], [1, 2])
    assert list(out.columns) == ["L.x", "L2.x", "L.y", "L2.y"]
    assert isinstance(panel.lag("x", [1]), pd.DataFrame)
    assert isinstance(panel.lag("x", 1), pd.Series)
    assert list(panel.operate("L(1/2).x", "L2.x", "x").columns) == ["L.x", "L2.x", "x"]


def test_unknown_and_non_numeric_variables_are_refused() -> None:
    frame = pd.DataFrame(
        {"id": [1, 1], "t": [1, 2], "x": [1.0, 2], "s": ["a", "b"], "b": [True, False]}
    )
    panel = xtset(frame, "id", "t")
    with pytest.raises(DataError, match="not numeric") as info:
        panel.lag("s")
    assert info.value.key == "panel.not_numeric"
    with pytest.raises(DataError):
        panel.lag("nope")
    same(panel.lag("b"), [NAN, 1.0])  # booleans are numbers


def test_integer_and_nullable_columns() -> None:
    frame = pd.DataFrame(
        {"id": [1, 1, 1], "t": [1, 2, 3], "k": pd.array([1, None, 3], dtype="Int64")}
    )
    same(xtset(frame, "id", "t").lag("k"), [NAN, 1.0, NAN])


def test_a_selection_keeps_the_calendar() -> None:
    panel = make([1, 2, 3, 4], [1.0, 2, 3, 4]).select("x != 3")
    same(panel.lag("x"), [NAN, 1.0, NAN, NAN])  # the removed row is missing; 4 has no lag


# ---------------------------------------------------------------- registry


def test_the_operators_are_registered_with_their_provenance() -> None:
    assert {"lag", "lead", "difference", "seasonal"} <= set(OPERATORS.names())
    assert OPERATORS.spec("lag").stata == "L.x, L2.x, ..."
