from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxtpanel._config import option_context
from pyxtpanel._errors import DataError, DeclarationError, SpecificationError, TimeAxisError
from pyxtpanel.panel.base import PanelBase

INVEST = Path(__file__).resolve().parents[2] / "validation" / "data" / "invest.dta"


def small() -> pd.DataFrame:
    """Unit 1 complete (t=1..4), unit 2 skips t=3, unit 3 enters late; rows shuffled."""
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 1, 2, 2, 2, 3, 3],
            "t": [1, 2, 3, 4, 1, 2, 4, 2, 3],
            "x": [1.0, 2, 3, 4, 5, 5, 7, 2, 2],
        }
    )
    return frame.sample(frac=1, random_state=3).reset_index(drop=True)


# ---------------------------------------------------------------- declaration


def test_declaration_builds_a_sorted_index_of_a_shuffled_frame() -> None:
    frame = small()
    panel = PanelBase.from_frame(frame, "id", "t")
    assert (panel.n_units, panel.n_obs, len(panel)) == (3, 9, 9)
    assert panel.index.offsets.tolist() == [0, 4, 7, 9]
    assert panel.index.unit.tolist() == [0, 0, 0, 0, 1, 1, 1, 2, 2]
    assert panel.index.rank.tolist() == [0, 1, 2, 3, 0, 1, 3, 1, 2]
    # source_row maps every sorted row back to the user's frame
    back = frame.iloc[panel.index.source_row]
    assert back["id"].tolist() == [1, 1, 1, 1, 2, 2, 2, 3, 3]
    assert back["t"].tolist() == [1, 2, 3, 4, 1, 2, 4, 2, 3]
    assert panel.counts.tolist() == [4, 3, 2] and panel.unit_labels.tolist() == [1, 2, 3]


def test_structure_properties_match_stata_xtset() -> None:
    panel = PanelBase.from_frame(small(), "id", "t")
    assert panel.balanced == "unbalanced" and panel.has_gaps and panel.n_gaps == 1
    assert (panel.t_min, panel.t_max, panel.delta, panel.n_periods) == (1, 4, 1, 4)
    assert panel.t_mean == 3.0
    text = repr(panel)
    assert "Panel variable: id (unbalanced)" in text
    assert "Time variable: t, 1 to 4, but with a gap" in text
    assert "Delta: 1 unit" in text


def test_the_users_frame_is_never_modified() -> None:
    frame = small()
    before = frame.copy(deep=True)
    panel = PanelBase.from_frame(frame, "id", "t")
    panel.to_frame(order="sorted", index="keys")
    pd.testing.assert_frame_equal(frame, before)
    assert panel.frame is not frame


@pytest.mark.parametrize(
    "ids",
    [
        ["b", "b", "a", "a"],
        pd.Categorical(["b", "b", "a", "a"]),
        pd.array([2, 2, 1, 1], dtype="Int64"),
        [2.0, 2.0, 1.0, 1.0],
    ],
)
def test_identifier_types(ids) -> None:  # type: ignore[no-untyped-def]
    frame = pd.DataFrame({"i": ids, "t": [1, 2, 1, 2], "x": [1, 2, 3, 4]})
    panel = PanelBase.from_frame(frame, "i", "t")
    assert panel.n_units == 2 and panel.index.offsets.tolist() == [0, 2, 4]
    assert panel.index.source_row.tolist() == [2, 3, 0, 1]  # the unit sorted first is the second


def test_mixed_integer_and_text_identifiers_sort_numbers_first() -> None:
    frame = pd.DataFrame({"i": [1, "a", 1, "a"], "t": [1, 1, 2, 2]})
    panel = PanelBase.from_frame(frame, "i", "t")
    assert panel.unit_labels.tolist() == [1, "a"]


def test_identifiers_that_cannot_be_ordered_are_reported() -> None:
    frame = pd.DataFrame({"i": [(1, 2), {"a": 1}, (1, 2), {"a": 1}], "t": [1, 1, 2, 2]})
    with pytest.raises((DeclarationError, TypeError)):
        PanelBase.from_frame(frame, "i", "t")


def test_two_level_index_is_read_automatically() -> None:
    frame = small().set_index(["id", "t"])
    panel = PanelBase.from_frame(frame)
    assert (panel.id_name, panel.time_name, panel.n_units) == ("id", "t", 3)
    unnamed = frame.rename_axis([None, None])
    assert PanelBase.from_frame(unnamed).n_obs == 9
    back = PanelBase.from_frame(unnamed).to_frame(order="sorted", index="keys")
    assert back.index.names == ["id", "time"]


def test_keys_can_be_index_levels_by_name() -> None:
    frame = small().set_index("id")
    panel = PanelBase.from_frame(frame, "id", "t")
    assert panel.n_units == 3
    out = panel.to_frame(order="sorted", index="keys")
    assert out.index.names == ["id", "t"] and out.index[0] == (1, 1)


def test_one_key_alone_or_a_flat_index_is_a_specification_error() -> None:
    frame = small()
    with pytest.raises(SpecificationError, match="both id= and time=") as info:
        PanelBase.from_frame(frame, "id")
    assert info.value.stata_code == 198
    with pytest.raises(SpecificationError, match="two levels"):
        PanelBase.from_frame(frame)


def test_unknown_variable_suggests_close_names() -> None:
    with pytest.raises(DataError, match=r"'idd' not found.*Did you mean 'id'.*Available") as info:
        PanelBase.from_frame(small(), "idd", "t")
    assert info.value.key == "panel.variable_not_found"


def test_duplicated_key_column_is_refused() -> None:
    frame = pd.DataFrame([[1, 1, 5], [1, 2, 6]], columns=["id", "t", "t"])
    with pytest.raises(DeclarationError, match="more than once"):
        PanelBase.from_frame(frame, "id", "t")


# ---------------------------------------------------------------- duplicates and missing keys


def test_repeated_pairs_are_refused_like_stata_451() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1, 2, 2], "t": [1, 2, 2, 1, 1], "x": range(5)})
    with pytest.raises(DeclarationError, match=r"2 repeated pair.*\(id=1, t=2\)") as info:
        PanelBase.from_frame(frame, "id", "t")
    assert info.value.stata_code == 451 and info.value.key == "panel.duplicates"


def test_missing_keys_are_refused_with_the_rows() -> None:
    frame = pd.DataFrame({"id": [1, 1, None, 2], "t": [1, 2, 1, 1]})
    with pytest.raises(DeclarationError, match=r"'id' has 1 missing.*rows: \[2\]") as info:
        PanelBase.from_frame(frame, "id", "t")
    assert info.value.stata_code == 459
    frame = pd.DataFrame({"id": [1, 1, 2], "t": [1, np.nan, 1]})
    with pytest.raises(DeclarationError, match="'t' has 1 missing"):
        PanelBase.from_frame(frame, "id", "t")


def test_missing_keys_drop_reproduces_stata_and_counts() -> None:
    frame = pd.DataFrame({"id": [1, 1, None, 2, 2], "t": [1, 2, 1, np.nan, 2], "x": range(5)})
    panel = PanelBase.from_frame(frame, "id", "t", missing_keys="drop")
    assert (panel.n_obs, panel.n_dropped, panel.n_units) == (3, 2, 2)
    assert panel.index.source_row.tolist() == [0, 1, 4]
    assert "Ignored: 2 observation(s)" in repr(panel)
    assert panel.to_frame()["x"].tolist() == [0, 1, 4]


def test_all_rows_unusable_is_error_2000() -> None:
    frame = pd.DataFrame({"id": [None, None], "t": [1, 2]})
    with pytest.raises(DeclarationError) as info:
        PanelBase.from_frame(frame, "id", "t", missing_keys="drop")
    assert info.value.stata_code == 2000
    with pytest.raises(DeclarationError):
        PanelBase.from_frame(pd.DataFrame({"id": [], "t": []}), "id", "t")


def test_invalid_missing_keys_option() -> None:
    with pytest.raises(SpecificationError, match="missing_keys"):
        PanelBase.from_frame(small(), "id", "t", missing_keys="ignore")  # type: ignore[arg-type]


def test_time_error_names_the_frame_row_not_the_position_after_dropping() -> None:
    frame = pd.DataFrame({"id": [None, 1, 1], "t": [1.0, 1.0, 2.5]})
    with pytest.raises(TimeAxisError, match=r"row 2 \(id 1"):
        PanelBase.from_frame(frame, "id", "t", missing_keys="drop")


def test_messages_follow_the_language_option() -> None:
    frame = pd.DataFrame({"id": [1, 1], "t": [1, 1]})
    with option_context(lang="fr"), pytest.raises(DeclarationError, match="Observations répétées"):
        PanelBase.from_frame(frame, "id", "t")


# ---------------------------------------------------------------- balance


def balanced_frame(periods_by_unit: dict[int, list[int]]) -> pd.DataFrame:
    rows = [(u, t) for u, ts in periods_by_unit.items() for t in ts]
    return pd.DataFrame(rows, columns=["id", "t"])


def test_strongly_weakly_and_unbalanced() -> None:
    strong = PanelBase.from_frame(balanced_frame({1: [1, 2, 3], 2: [1, 2, 3]}), "id", "t")
    weak = PanelBase.from_frame(balanced_frame({1: [1, 2], 2: [2, 3]}), "id", "t")
    unbalanced = PanelBase.from_frame(balanced_frame({1: [1, 2, 3], 2: [1, 2]}), "id", "t")
    assert strong.balanced == "strongly balanced" and not strong.has_gaps
    assert weak.balanced == "weakly balanced"
    assert unbalanced.balanced == "unbalanced"


def test_balance_summary_names_the_first_incomplete_unit() -> None:
    panel = PanelBase.from_frame(balanced_frame({1: [1, 2, 3, 4], 2: [1, 2], 3: [1]}), "id", "t")
    summary = panel.balance_summary
    assert (summary.n_units, summary.n_periods, summary.n_incomplete) == (3, 4, 2)
    assert summary.first_incomplete == 1 and summary.missing_ranks == (2, 3)


def test_late_entry_is_not_a_gap() -> None:
    panel = PanelBase.from_frame(balanced_frame({1: [1, 2, 3], 2: [3]}), "id", "t")
    assert not panel.has_gaps and panel.balanced == "unbalanced"


# ---------------------------------------------------------------- time variants


def test_dates_with_a_frequency() -> None:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 2],
            "d": pd.to_datetime(["2020-01-01", "2020-04-01", "2020-10-01", "2020-01-01"]),
        }
    )
    panel = PanelBase.from_frame(frame, "id", "d")
    assert panel.axis.freq == "Q" and panel.has_gaps
    assert repr(panel).count("2020q1") == 1 and "2020q4" in repr(panel)
    assert panel.t_min == pd.Timestamp("2020-01-01")


def test_delta_and_off_grid_time() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1], "t": [1, 3, 7]})
    panel = PanelBase.from_frame(frame, "id", "t", delta=2)
    assert panel.has_gaps and panel.n_gaps == 1 and "Delta: 2 units" in repr(panel)
    with pytest.raises(TimeAxisError, match="not on the calendar grid"):
        PanelBase.from_frame(pd.DataFrame({"id": [1, 1], "t": [1, 2]}), "id", "t", delta=2)


# ---------------------------------------------------------------- back to a frame


def test_to_frame_orders_and_index() -> None:
    frame = small()
    panel = PanelBase.from_frame(frame, "id", "t")
    original = panel.to_frame()
    pd.testing.assert_frame_equal(original, frame)
    ordered = panel.to_frame(order="sorted")
    assert ordered["id"].tolist() == [1, 1, 1, 1, 2, 2, 2, 3, 3]
    keyed = panel.to_frame(order="sorted", index="keys")
    assert keyed.index.names == ["id", "t"] and list(keyed.columns) == ["x"]
    assert keyed.index[4] == (2, 1)


# ---------------------------------------------------------------- Hansen panel (against Stata)


@pytest.mark.skipif(not INVEST.exists(), reason="validation data not present")
def test_hansen_panel_matches_stata_xtset_and_xtdescribe() -> None:
    frame = pd.read_stata(INVEST)
    panel = PanelBase.from_frame(frame, "firm", "year")
    # Stata: xtset firm year -> strongly balanced, 1973 to 1987, Delta 1 unit; xtdescribe: n = 565
    assert panel.balanced == "strongly balanced" and not panel.has_gaps
    assert (panel.n_units, panel.n_obs, panel.t_min, panel.t_max) == (565, 8475, 1973, 1987)
    assert (panel.counts == 15).all() and panel.axis.delta_text == "1 unit"


# ---------------------------------------------------------------- properties


@pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")
@settings(max_examples=60, deadline=None)
@given(
    data=st.lists(
        st.tuples(st.integers(0, 5), st.integers(0, 9)), min_size=1, max_size=40, unique=True
    ),
    seed=st.integers(0, 1000),
)
def test_row_order_does_not_change_the_sorted_index(data, seed) -> None:  # type: ignore[no-untyped-def]
    frame = pd.DataFrame(data, columns=["id", "t"])
    shuffled = frame.sample(frac=1, random_state=seed).reset_index(drop=True)
    a = PanelBase.from_frame(frame, "id", "t")
    b = PanelBase.from_frame(shuffled, "id", "t")
    assert (a.index.offsets == b.index.offsets).all()
    assert (a.index.unit == b.index.unit).all() and (a.index.rank == b.index.rank).all()
    assert a.balanced == b.balanced and a.n_gaps == b.n_gaps
    # and the unit labels are renamed consistently (an offset of 100 keeps the order)
    c = PanelBase.from_frame(frame.assign(id=frame["id"] + 100), "id", "t")
    assert (a.index.offsets == c.index.offsets).all() and (a.index.rank == c.index.rank).all()


@settings(max_examples=80, deadline=None)
@given(
    data=st.lists(
        st.tuples(st.integers(0, 8), st.integers(0, 15)), min_size=1, max_size=60, unique=True
    ),
    seed=st.integers(0, 1000),
)
def test_the_order_found_by_the_duplicate_check_is_the_one_of_a_separate_sort(data, seed) -> None:  # type: ignore[no-untyped-def]
    """Declaring sorts once (I.2.1): the shared order must equal ``lexsort`` on (unit, rank)."""
    from pyxtpanel.panel.checks import check_unique_pairs
    from pyxtpanel.panel.index import build_index

    pairs = np.array(data)[np.random.default_rng(seed).permutation(len(data))]
    unit, rank = pairs[:, 0].astype(np.int64), pairs[:, 1].astype(np.int64)
    shown = (unit, rank)
    order = check_unique_pairs(unit, rank, "id", "t", shown)
    assert (order == np.lexsort((rank, unit))).all()
    labels = np.arange(unit.max() + 1)
    positions = np.arange(unit.size, dtype=np.int64)
    shared = build_index(unit, rank, positions, labels, order)
    alone = build_index(unit, rank, positions, labels)
    for name in ("unit_labels", "offsets", "unit", "rank", "source_row"):
        assert (getattr(shared, name) == getattr(alone, name)).all()


def test_alignment_of_float_results_puts_missing_where_rows_were_set_aside() -> None:
    frame = pd.DataFrame({"id": [1, 1, 2, 2, np.nan], "t": [1, 2, 1, 2, 1], "x": [1.0, 2, 3, 4, 5]})
    panel = PanelBase.from_frame(frame, "id", "t", missing_keys="drop")
    out = panel.align(np.array([10.0, 20.0, 30.0, 40.0]), "out")
    assert out.tolist()[:4] == [10.0, 20.0, 30.0, 40.0] and np.isnan(out.iloc[4])
    assert out.name == "out" and out.index.equals(frame.index)
    ints = panel.align(np.array([1, 2, 3, 4]), "k")  # integers keep the reindex path (float, NaN)
    assert ints.iloc[:4].tolist() == [1.0, 2.0, 3.0, 4.0] and np.isnan(ints.iloc[4])
