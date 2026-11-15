import numpy as np
import pandas as pd
import pytest

from pyxtpanel._errors import DataError, DeclarationError, SpecificationError, TimeAxisError
from pyxtpanel.panel.base import PanelBase
from pyxtpanel.panel.gaps import find_duplicates


def probe() -> pd.DataFrame:
    """The panel used in the Stata probe (validation/stata/probe_xt.do), in shuffled order."""
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 1, 2, 2, 2, 3, 3],
            "t": [1, 2, 3, 4, 1, 2, 4, 2, 3],
            "x": [1.0, 2, 3, 4, 5, 5, 7, 2, 2],
            "z": [10.0, 11, 12, 13, np.nan, 15, 17, 20, 21],
        }
    )
    return frame.sample(frac=1, random_state=11).reset_index(drop=True)


# ---------------------------------------------------------------- selection


def test_where_mask_and_expression() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    by_expr = panel.select("x > 2")
    assert by_expr.n_obs == 5 and by_expr.n_units == 2  # unit 3 only has x=2
    mask = (probe()["id"] == 2).to_numpy()
    assert panel.select(mask).n_units == 1
    assert panel.n_obs == 9  # the original panel is untouched


def test_missing_in_condition_counts_as_false() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    assert panel.select("z > 0").n_obs == 8  # the row with missing z is out


def test_where_must_be_boolean_and_full_length() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    with pytest.raises(DataError, match="one True/False per row"):
        panel.select(np.array([True, False]))
    with pytest.raises(DataError, match="one True/False per row"):
        panel.select("x")


def test_units_periods_and_complete() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    assert panel.select(units=[1, 3]).unit_labels.tolist() == [1, 3]
    assert panel.select(periods=[1, 2]).n_obs == 5
    assert panel.select(periods=slice(2, 3)).n_obs == 5
    assert panel.select(periods=slice(3, None)).counts.tolist() == [2, 1, 1]
    assert panel.select(complete=["z"]).n_obs == 8
    assert panel.select(complete=["x", "z"], units=[2]).n_obs == 2


def test_selection_keeps_declaration_and_updates_summary() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    sub = panel.select(units=[1, 3])
    assert (sub.id_name, sub.time_name, sub.axis) == (panel.id_name, panel.time_name, panel.axis)
    assert sub.balanced == "unbalanced" and (sub.t_min, sub.t_max) == (1, 4)
    late = panel.select(periods=slice(2, 3))
    assert (late.t_min, late.t_max) == (2, 3) and "2 to 3" in repr(late)
    only_one = panel.select(units=[1]).select(periods=slice(1, 3))
    assert only_one.balanced == "strongly balanced"


def test_selecting_nothing_or_unknown_inputs() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    with pytest.raises(DeclarationError) as info:
        panel.select(units=[99])
    assert info.value.stata_code == 2000
    with pytest.raises(DataError, match="'zz' not found"):
        panel.select(complete=["zz"])
    with pytest.raises(TimeAxisError):
        PanelBase.from_frame(probe(), "id", "t", delta=1).select(periods=[1.5])


def test_infinite_values_are_not_complete() -> None:
    frame = probe()
    frame.loc[0, "x"] = np.inf
    assert PanelBase.from_frame(frame, "id", "t").select(complete=["x"]).n_obs == 8


# ---------------------------------------------------------------- gaps and duplicates


def test_gaps_table() -> None:
    table = PanelBase.from_frame(probe(), "id", "t").gaps()
    assert table.index.name == "id"
    assert table["first"].tolist() == [1, 1, 2] and table["last"].tolist() == [4, 4, 3]
    assert table["n_obs"].tolist() == [4, 3, 2]
    assert table["n_breaks"].tolist() == [0, 1, 0] and table["n_missing"].tolist() == [0, 1, 0]
    only = PanelBase.from_frame(probe(), "id", "t").gaps(only_with_gaps=True)
    assert only.index.tolist() == [2]


def test_gaps_counts_several_breaks_and_missing_periods() -> None:
    frame = pd.DataFrame({"id": [1] * 4, "t": [1, 3, 4, 8]})
    row = PanelBase.from_frame(frame, "id", "t").gaps().iloc[0]
    assert (row["n_breaks"], row["n_missing"], row["n_obs"]) == (2, 4, 4)


def test_find_duplicates_lists_all_repeated_rows() -> None:
    frame = pd.DataFrame(
        {"id": [1, 1, 2, 2, 2, 3, None], "t": [1, 1, 5, 5, 5, 1, 1], "x": range(7)},
        index=list("abcdefg"),
    )
    out = find_duplicates(frame, "id", "t")
    assert out.index.tolist() == ["a", "b", "c", "d", "e"]
    assert out["n_repeats"].tolist() == [2, 2, 3, 3, 3]
    assert find_duplicates(frame.dropna().drop_duplicates(["id", "t"]), "id", "t").empty


# ---------------------------------------------------------------- tsfill


def test_tsfill_matches_stata() -> None:
    # Stata (probe_xt.do): tsfill adds (2,3); tsfill, full also adds (3,1) and (3,4)
    panel = PanelBase.from_frame(probe(), "id", "t")
    filled = panel.tsfill(indicator="_new")
    data = filled.to_frame(order="sorted")
    assert data[["id", "t"]].values.tolist() == [
        [1, 1], [1, 2], [1, 3], [1, 4], [2, 1], [2, 2], [2, 3], [2, 4], [3, 2], [3, 3]
    ]  # fmt: skip
    added = data[data["_new"]]
    assert added[["id", "t"]].values.tolist() == [[2, 3]]
    assert added[["x", "z"]].isna().all().all()
    assert not filled.has_gaps and filled.n_obs == 10

    full = panel.tsfill(full=True).to_frame(order="sorted")
    assert full[["id", "t"]].values.tolist() == ([[u, t] for u in (1, 2, 3) for t in (1, 2, 3, 4)])
    assert full.loc[(full["id"] == 3) & (full["t"] == 4), ["x", "z"]].isna().all().all()
    assert (
        PanelBase.from_frame(probe(), "id", "t").tsfill(full=True).balanced == "strongly balanced"
    )


def test_tsfill_with_delta_dates_and_index_keys() -> None:
    dates = pd.to_datetime(["2020-01-01", "2020-07-01", "2021-07-01"])
    frame = pd.DataFrame({"id": [1, 1, 1], "d": dates, "x": [1.0, 2, 3]}).set_index("id")
    panel = PanelBase.from_frame(frame, "id", "d", freq="H")
    filled = panel.tsfill()
    out = filled.to_frame(order="sorted")
    assert len(out) == 4 and out["d"].tolist()[2] == pd.Timestamp("2021-01-01")
    delta = PanelBase.from_frame(pd.DataFrame({"id": [1, 1], "t": [1, 7]}), "id", "t", delta=3)
    assert delta.tsfill().to_frame(order="sorted")["t"].tolist() == [1, 4, 7]


def test_tsfill_indicator_collision_and_noop() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    with pytest.raises(SpecificationError, match="already exists"):
        panel.tsfill(indicator="x")
    complete = PanelBase.from_frame(pd.DataFrame({"id": [1, 1], "t": [1, 2]}), "id", "t")
    assert complete.tsfill(indicator="_f").to_frame()["_f"].tolist() == [False, False]


# ---------------------------------------------------------------- balance


def test_balance_drop_and_min_obs() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    assert panel.balance().unit_labels.tolist() == [1]
    assert panel.balance().balanced == "strongly balanced"
    assert panel.balance(min_obs=3).unit_labels.tolist() == [1]  # drop still needs all periods


def test_balance_largest_block() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    best = panel.balance("largest_block")
    # candidate blocks: unit 1 alone on 1..4 = 4 obs; units 1,2,3 on t=2..? unit 2 lacks 3;
    # units 1,2 on t=1..2 = 4 obs; units 1,2,3 on t=2 = 3 obs; units 1,3 on t=2..3 = 4 obs
    assert best.n_obs == 4 and best.balanced == "strongly balanced"
    assert best.unit_labels.tolist() == [1]  # ties go to the longer window


def test_balance_largest_block_prefers_more_observations() -> None:
    frame = pd.DataFrame(
        {
            "id": [1] * 5 + [2] * 3 + [3] * 3 + [4] * 3,
            "t": [*range(1, 6), 1, 2, 3, 2, 3, 4, 3, 4, 5],
        }
    )
    best = PanelBase.from_frame(frame, "id", "t").balance("largest_block")
    assert best.n_obs == 6 and best.t_min >= 1  # e.g. 3 units x 2 periods or 2 units x 3 periods
    assert best.balanced == "strongly balanced"


def test_balance_errors() -> None:
    panel = PanelBase.from_frame(probe(), "id", "t")
    with pytest.raises(SpecificationError, match="method must be"):
        panel.balance("fill")  # type: ignore[arg-type]
    with pytest.raises(DeclarationError):
        panel.balance(min_obs=10)


# ---------------------------------------------------------------- pandas bridges and alignment


def test_multiindex_roundtrip_and_matches() -> None:
    frame = probe()
    panel = PanelBase.from_frame(frame, "id", "t")
    wide = panel.to_multiindex()
    assert wide.index.names == ["id", "t"] and wide.index[0] == (1, 1)
    again = PanelBase.from_multiindex(wide)
    assert again.n_obs == panel.n_obs and again.balanced == panel.balanced
    assert panel.matches(frame)
    edited = frame.copy()
    edited.loc[0, "x"] = -1.0
    assert not panel.matches(edited)


def test_later_edits_of_the_original_do_not_reach_the_panel() -> None:
    frame = probe()
    panel = PanelBase.from_frame(frame, "id", "t")
    frame.loc[0, "x"] = -999.0
    assert (panel.frame["x"] != -999.0).all() and not panel.matches(frame)


def test_align_puts_values_back_on_the_original_rows() -> None:
    frame = probe()
    frame.index = list("abcdefghi")
    panel = PanelBase.from_frame(frame, "id", "t")
    sorted_x = panel.frame["x"].to_numpy()[panel.index.source_row]
    back = panel.align(sorted_x, "x_again")
    pd.testing.assert_series_equal(back, frame["x"].rename("x_again"))
    two = panel.align(np.column_stack([sorted_x, sorted_x * 2]))
    assert two.shape == (9, 2) and (two[1] == frame["x"] * 2).all()


def test_align_leaves_rows_outside_the_panel_missing() -> None:
    frame = probe()
    panel = PanelBase.from_frame(frame, "id", "t").select(units=[1])
    out = panel.align(np.arange(panel.n_obs, dtype=float))
    assert out.notna().sum() == 4 and out[frame["id"] != 1].isna().all()
    with pytest.raises(ValueError, match="expected 4 values"):
        panel.align(np.arange(3.0))


def test_index_with_repeated_labels_is_supported() -> None:
    frame = probe()
    frame.index = [0] * len(frame)
    panel = PanelBase.from_frame(frame, "id", "t")
    assert panel.n_obs == 9
    back = panel.align(np.arange(9.0))
    assert len(back) == 9 and sorted(back.tolist()) == list(map(float, range(9)))


def test_other_dtypes_in_the_frame() -> None:
    frame = pd.DataFrame(
        {
            "id": ["a", "a", "b"],
            "t": [1, 2, 1],
            "s": pd.array(["u", "v", None], dtype="string"),
            "c": pd.Categorical(["p", "q", "p"]),
            "n": pd.array([1, None, 3], dtype="Int64"),
        }
    )
    panel = PanelBase.from_frame(frame, "id", "t")
    assert panel.select(complete=["s"]).n_obs == 2
    assert panel.select(complete=["n"]).n_obs == 2
    assert panel.select(complete=["c"]).n_obs == 3


def test_linearmodels_bridge_by_duck_typing() -> None:
    class Fake:
        dataframe = probe().set_index(["id", "t"])

    panel = PanelBase.from_linearmodels(Fake())
    assert panel.n_units == 3 and panel.n_obs == 9
    assert PanelBase.from_linearmodels(probe().set_index(["id", "t"])).n_obs == 9
    with pytest.raises(SpecificationError, match="indexed by"):
        PanelBase.from_linearmodels(probe())
    pytest.importorskip("linearmodels")
    assert panel.to_linearmodels() is not None


def test_with_frame_reuses_the_sorted_index_for_new_values() -> None:
    frame = probe()
    panel = PanelBase.from_frame(frame, "id", "t")
    changed = frame.assign(x=frame["x"] * 10)
    again = panel.with_frame(changed)
    assert again.index is panel.index and again.axis is panel.axis
    assert again.frame["x"].tolist() == (frame["x"] * 10).tolist()
    assert panel.frame["x"].tolist() == frame["x"].tolist()  # the original panel is untouched


def test_with_frame_refuses_other_rows_or_keys() -> None:
    frame = probe()
    panel = PanelBase.from_frame(frame, "id", "t")
    with pytest.raises(DataError, match="8 rows instead of 9"):
        panel.with_frame(frame.iloc[:8])
    with pytest.raises(DataError, match="values differ") as info:
        panel.with_frame(frame.assign(t=frame["t"] + 1))
    assert info.value.key == "panel.with_frame_mismatch"
    indexed = frame.set_index(["id", "t"]).rename_axis([None, None])
    unnamed = PanelBase.from_frame(indexed)
    assert unnamed.with_frame(indexed.assign(x=0.0)).n_obs == 9
    with pytest.raises(DataError):
        unnamed.with_frame(indexed.iloc[::-1])
