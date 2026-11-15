"""Wide, long and array forms and duplicates (K.5): by hand, round trips, loop oracle."""

from __future__ import annotations

import math
import statistics

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import pyxtpanel as xt
from pyxtpanel._errors import SpecificationError

quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")
nan = math.nan


def hand() -> xt.PanelData:
    rows = [
        ("A", 1990, 1.0, 10.0),
        ("A", 1991, 2.0, 20.0),
        ("A", 1992, 3.0, 30.0),
        ("B", 1990, 4.0, 40.0),
        ("B", 1992, 6.0, 60.0),
    ]
    return xt.xtset(pd.DataFrame(rows, columns=["id", "t", "x", "z"]), "id", "t")


# ----------------------------------------------------------------- to_wide / to_long


def test_to_wide_by_hand() -> None:
    wide = hand().to_wide(["x", "z"])
    assert wide.columns.tolist() == ["x1990", "x1991", "x1992", "z1990", "z1991", "z1992"]
    assert wide.index.tolist() == ["A", "B"] and wide.index.name == "id"
    np.testing.assert_array_equal(
        wide.to_numpy(),
        [[1, 2, 3, 10, 20, 30], [4, nan, 6, 40, nan, 60]],
    )
    assert hand().to_wide(["x"], sep="_").columns.tolist() == ["x_1990", "x_1991", "x_1992"]


def test_to_wide_defaults_to_the_numeric_variables_and_labels_dates() -> None:
    assert hand().to_wide().columns.tolist()[:2] == ["x1990", "x1991"]
    quarters = pd.date_range("2020-01-01", periods=3, freq="QS")
    frame = pd.DataFrame({"id": 1, "t": quarters, "x": [1.0, 2, 3]})
    wide = xt.xtset(frame, "id", "t", freq="Q").to_wide(["x"])
    assert wide.columns.tolist() == ["x2020q1", "x2020q2", "x2020q3"]


def test_to_long_by_hand() -> None:
    wide = pd.DataFrame(
        {
            "id": ["A", "B"],
            "x1990": [1.0, 4.0],
            "x1991": [2.0, nan],
            "z1990": [10.0, 40.0],
            "g": ["u", "v"],
        }
    )
    long = xt.to_long(wide, ["x", "z"], id="id", time="year")
    assert long.columns.tolist() == ["id", "year", "x", "z", "g"]
    assert long["id"].tolist() == ["A", "A", "B", "B"]
    assert long["year"].tolist() == [1990, 1991, 1990, 1991]
    np.testing.assert_array_equal(long["x"].to_numpy(), [1, 2, 4, nan])
    np.testing.assert_array_equal(long["z"].to_numpy(), [10, nan, 40, nan])
    assert long["g"].tolist() == ["u", "u", "v", "v"]
    kept = xt.to_long(wide, ["x", "z"], id="id", time="year", drop_empty=True)
    assert len(kept) == 3  # B in 1991 has nothing observed


def test_to_long_from_the_index_and_with_named_suffixes() -> None:
    wide = pd.DataFrame(
        {"x_q1": [1.0, 2.0], "x_q2": [3.0, 4.0]}, index=pd.Index([7, 8], name="firm")
    )
    long = xt.to_long(wide, ["x"], sep="_", times={"q1": 1, "q2": 2})
    assert long.columns.tolist() == ["firm", "t", "x"]
    assert long["firm"].tolist() == [7, 7, 8, 8] and long["x"].tolist() == [1.0, 3.0, 2.0, 4.0]
    with pytest.raises(SpecificationError) as bad:
        xt.to_long(wide, ["x"], sep="_")  # the suffixes are not integers
    assert bad.value.stata_code == 198
    with pytest.raises(SpecificationError):
        xt.to_long(wide, ["x"], id="missing", sep="_")


@quiet
@settings(max_examples=60, deadline=None)
@given(
    data=st.dictionaries(
        st.integers(0, 6),
        st.dictionaries(st.integers(1990, 2000), st.integers(-9, 9).map(float), min_size=1),
        min_size=1,
    )
)
def test_wide_then_long_gives_back_the_panel(data) -> None:  # type: ignore[no-untyped-def]
    rows = [(u, t, x) for u, rs in data.items() for t, x in rs.items()]
    frame = pd.DataFrame(rows, columns=["id", "t", "x"])
    wide = xt.xtset(frame, "id", "t").to_wide(["x"]).reset_index()
    back = xt.to_long(wide, ["x"], id="id", time="t", drop_empty=True)
    expected = frame.sort_values(["id", "t"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(back, expected, check_dtype=False)


# ----------------------------------------------------------------- to_array3d


def test_to_array3d_by_hand() -> None:
    cube = hand().to_array3d(["x", "z"])
    assert cube.values.shape == (2, 3, 2)  # variable, period, unit
    assert cube.variables == ("x", "z") and cube.units.tolist() == ["A", "B"]
    assert cube.periods.tolist() == [1990, 1991, 1992]
    np.testing.assert_array_equal(cube.values[0], [[1, 4], [2, nan], [3, 6]])
    np.testing.assert_array_equal(cube.observed, [[True, True], [True, False], [True, True]])
    assert np.isnan(cube.values[1, 1, 1])


def test_the_cube_covers_periods_nobody_is_observed_in_and_agrees_with_the_contract() -> None:
    frame = pd.DataFrame({"id": [1, 1, 2, 2], "t": [1, 4, 1, 4], "x": [1.0, 2, 3, 4]})
    cube = xt.xtset(frame, "id", "t").to_array3d(["x"])
    assert cube.periods.tolist() == [1, 2, 3, 4]
    assert not cube.observed[1:3].any()
    balanced = xt.xtset(frame.assign(t=[1, 2, 1, 2]), "id", "t")
    contract = balanced.arrays(["x"]).as_cube()  # (units, periods, variables)
    np.testing.assert_array_equal(balanced.to_array3d(["x"]).values[0].T, contract[:, :, 0])


# ----------------------------------------------------------------- duplicates


def dupes() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": [1, 1, 1, 2, 2, 1, 3.0],
            "t": [1, 1, 2, 1, 1, 1, np.nan],
            "x": [1.0, 3.0, 5.0, nan, 7.0, 5.0, 9.0],
            "s": ["a", "b", "c", "d", "e", "f", "g"],
        }
    )


def test_first_and_last_keep_one_row_whole() -> None:
    first = xt.resolve_duplicates(dupes(), "id", "t", how="first")
    assert first.frame["s"].tolist() == ["a", "c", "d", "g"]
    assert first.n_rows_in == 7 and first.n_rows_out == 4
    last = xt.resolve_duplicates(dupes(), "id", "t", how="last")
    assert last.frame["s"].tolist() == ["c", "e", "f", "g"]
    assert first.groups.to_dict("records") == [
        {"id": 1.0, "t": 1.0, "n_rows": 3, "differing": "x, s"},
        {"id": 2.0, "t": 1.0, "n_rows": 2, "differing": "x, s"},
    ]


def test_statistics_replace_the_numeric_columns_and_ignore_missing_values() -> None:
    mean = xt.resolve_duplicates(dupes(), "id", "t", how="mean").frame
    assert mean["x"].tolist() == [3.0, 5.0, 7.0, 9.0] and mean["s"].tolist() == ["a", "c", "d", "g"]
    assert mean["id"].tolist() == [1, 1, 2, 3] and mean["t"].isna().tolist() == [
        False,
        False,
        False,
        True,
    ]
    assert xt.resolve_duplicates(dupes(), "id", "t", how="sum").frame["x"].tolist() == [9, 5, 7, 9]
    allmissing = pd.DataFrame({"id": [1, 1], "t": [1, 1], "x": [nan, nan]})
    assert xt.resolve_duplicates(allmissing, "id", "t", how="sum").frame["x"].isna().all()


def test_the_resolved_frame_can_be_declared_and_the_report_reads() -> None:
    result = xt.resolve_duplicates(dupes().iloc[:6], "id", "t", how="mean")
    assert xt.xtset(result.frame, "id", "t").n_obs == 3
    assert "2 repeated pairs, 3 rows removed" in result.summary("en")
    assert "2 paires répétées" in result.summary("fr") and "3 lignes supprimées" in result.summary(
        "fr"
    )
    assert xt.resolve_duplicates(hand().to_frame(), "id", "t").summary("en") == "No repeated pair."


def test_a_two_level_index_and_a_bad_rule() -> None:
    frame = dupes().iloc[:3].set_index(["id", "t"])
    result = xt.resolve_duplicates(frame, how="mean")
    assert len(result.frame) == 2
    with pytest.raises(SpecificationError) as bad:
        xt.resolve_duplicates(dupes(), "id", "t", how="mode")  # type: ignore[arg-type]
    assert bad.value.stata_code == 198


def loop_resolve(rows: list[tuple[int, int, float]], how: str) -> list[tuple[int, int, float]]:
    order: list[tuple[int, int]] = []
    seen: dict[tuple[int, int], list[float]] = {}
    for u, t, x in rows:
        if (u, t) not in seen:
            order.append((u, t))
            seen[(u, t)] = []
        seen[(u, t)].append(x)
    out = []
    for key in order:
        values = seen[key]
        good = [v for v in values if not math.isnan(v)]
        if how == "first":
            value = values[0]
        elif how == "last":
            value = values[-1]
        elif not good:
            value = nan
        else:
            value = {
                "mean": statistics.fmean,
                "sum": sum,
                "min": min,
                "max": max,
                "median": statistics.median,
            }[how](good)
        out.append((key[0], key[1], value))
    return out


@quiet
@settings(max_examples=120, deadline=None)
@given(
    rows=st.lists(
        st.tuples(
            st.integers(0, 3),
            st.integers(0, 3),
            st.one_of(st.integers(-5, 5).map(float), st.just(nan)),
        ),
        min_size=1,
        max_size=14,
    ),
    how=st.sampled_from(["first", "last", "mean", "sum", "min", "max", "median"]),
)
def test_resolution_equals_the_loop_implementation(rows, how) -> None:  # type: ignore[no-untyped-def]
    frame = pd.DataFrame(rows, columns=["id", "t", "x"])
    result = xt.resolve_duplicates(frame, "id", "t", how=how)
    expected = loop_resolve(rows, how)
    # first/last keep rows in their original position; the statistics keep the first row's position
    got = sorted(zip(result.frame["id"], result.frame["t"], result.frame["x"], strict=True))
    want = sorted(expected)
    assert len(got) == len(want)
    for (gu, gt, gx), (wu, wt, wx) in zip(got, want, strict=True):
        assert (gu, gt) == (wu, wt)
        np.testing.assert_allclose(gx, wx, rtol=1e-12, equal_nan=True)
    assert not result.frame.duplicated(["id", "t"]).any()
    assert result.n_rows_in - result.n_rows_out == int(
        result.groups["n_rows"].sum() - len(result.groups)
    )
