"""Contract test of the surface offered to companion packages (B.1.2).

``pyxtpanel.extend``, ``panel.calendar`` and ``panel.term_arrays`` are what ``pyxtviz`` (and any
other package built on the foundation) relies on instead of private modules.  The tests check the
proofs the design asks for: every re-export is the internal object, the calendar agrees with the
internal axis, the term reader agrees with ``arrays`` and ``operate``, and importing the surface
does not pull in matplotlib.
"""

from __future__ import annotations

import inspect
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel import _i18n, _registry, _signature, _types, extend
from pyxtpanel._errors import DataError, TimeAxisError
from pyxtpanel.contract import CONTRACT_VERSION, PanelArrays
from pyxtpanel.time.calendar import Calendar

# ---------------------------------------------------------------- pyxtpanel.extend


def test_each_reexport_is_the_internal_object() -> None:
    assert extend.Registry is _registry.Registry
    assert extend.register_messages is _i18n.register
    assert extend.message is _i18n.message
    assert extend.message_keys is _i18n.all_keys
    assert extend.has_language is _i18n.has_language
    assert extend.message_fields is _i18n.fields_of
    assert extend.call_reason is _signature.reason
    assert extend.call_arguments is _signature.accepted
    assert extend.LANGUAGES is _types.LANGUAGES
    assert extend.Lang is _types.Lang
    assert extend.FloatArray is _types.FloatArray
    assert extend.IntArray is _types.IntArray
    assert extend.BoolArray is _types.BoolArray


def test_the_surface_is_the_documented_one() -> None:
    assert set(extend.__all__) == {
        "EXTEND_VERSION",
        "ExtendedPanelLike",
        "Registry",
        "register_messages",
        "message",
        "message_keys",
        "has_language",
        "message_fields",
        "call_reason",
        "call_arguments",
        "Lang",
        "LANGUAGES",
        "FloatArray",
        "IntArray",
        "BoolArray",
    }
    assert extend.EXTEND_VERSION == 1
    assert CONTRACT_VERSION == 1  # the array contract did not change


def test_importing_extend_does_not_import_matplotlib() -> None:
    code = (
        "import sys, pyxtpanel.extend; "
        "raise SystemExit(1 if any(m.split('.')[0] == 'matplotlib' for m in sys.modules) else 0)"
    )
    assert subprocess.run([sys.executable, "-c", code], check=False).returncode == 0


def test_a_companion_registers_its_own_messages_and_reads_a_bad_call() -> None:
    extend.register_messages(
        "test_ext",
        {"test_ext.hello": {"en": "hello {who}", "fr": "bonjour {who}"}},
    )
    assert extend.message("test_ext.hello", "fr", who="monde") == "bonjour monde"
    assert "test_ext.hello" in extend.message_keys()
    assert extend.has_language("test_ext.hello", "fr")
    assert extend.message_fields("test_ext.hello", "en") == {"who"}

    def draw(panel: object, variable: str, *, bins: int = 10) -> None: ...

    signature = inspect.signature(draw)
    with pytest.raises(TypeError) as error:
        signature.bind(object())
    assert extend.call_reason(signature, error.value) == "missing a required argument: 'variable'"
    assert extend.call_arguments(signature) == "variable, bins=10"


# ---------------------------------------------------------------- panel.calendar


def _integer_panel() -> xt.PanelData:
    frame = pd.DataFrame({"id": [1, 1, 1, 2, 2], "t": [1990, 1991, 1993, 1990, 1992], "x": 1.0})
    return xt.xtset(frame, "id", "t")


def _date_panel() -> xt.PanelData:
    dates = pd.to_datetime(["2020-01-31", "2020-02-29", "2020-04-30"])
    frame = pd.DataFrame({"id": [1, 1, 1], "t": dates, "x": [1.0, 2.0, 3.0]})
    return xt.xtset(frame, "id", "t")


def _period_panel() -> xt.PanelData:
    periods = pd.PeriodIndex(["2019Q1", "2019Q2", "2019Q4"], freq="Q")
    frame = pd.DataFrame({"id": [1, 1, 1], "t": periods, "x": [1.0, 2.0, 3.0]})
    return xt.xtset(frame, "id", "t")


PANELS = {"integer": _integer_panel, "datetime": _date_panel, "period": _period_panel}


@pytest.mark.parametrize("name", PANELS)
def test_the_calendar_agrees_with_the_internal_axis(name: str) -> None:
    panel = PANELS[name]()
    calendar, axis = panel.calendar, panel.axis
    assert isinstance(calendar, Calendar)
    assert (calendar.kind, calendar.freq, calendar.delta, calendar.n_steps) == (
        axis.kind,
        axis.freq,
        axis.delta,
        axis.n_steps,
    )
    positions = panel.arrays(["x"]).time_index
    assert calendar.values(positions).equals(axis.values(positions))
    assert calendar.labels(positions) == axis.labels(positions)
    everything = np.arange(axis.n_steps)
    assert calendar.labels(everything) == axis.labels(everything)


@pytest.mark.parametrize("name", ["datetime", "period"])
def test_datetimes_of_a_dated_calendar(name: str) -> None:
    panel = PANELS[name]()
    positions = np.arange(panel.calendar.n_steps)
    assert panel.calendar.datetimes(positions).equals(panel.axis.datetimes(positions))


@pytest.mark.parametrize("name", PANELS)
def test_positions_of_values_is_the_inverse_of_values(name: str) -> None:
    panel = PANELS[name]()
    positions = np.arange(panel.calendar.n_steps)
    assert np.array_equal(panel.calendar.positions_of(panel.calendar.values(positions)), positions)


def test_the_positions_of_the_gaps_are_the_ones_of_the_contract() -> None:
    panel = _integer_panel()
    arrays = panel.arrays(["x"])
    assert arrays.time_index.tolist() == [0, 1, 3, 0, 2]
    assert panel.calendar.values(arrays.time_index).tolist() == [1990, 1991, 1993, 1990, 1992]


@pytest.mark.parametrize("lang", ["en", "fr"])
def test_datetimes_of_an_integer_calendar_is_a_bilingual_error(lang: str) -> None:
    xt.set_options(lang=lang)
    try:
        with pytest.raises(TimeAxisError) as error:
            _integer_panel().calendar.datetimes([0])
    finally:
        xt.set_options(lang="en")
    assert error.value.key == "time.no_dates"
    assert ("no dates" if lang == "en" else "pas de dates") in str(error.value)


def test_the_calendar_is_read_only() -> None:
    calendar = _integer_panel().calendar
    with pytest.raises(AttributeError):
        calendar.delta = 3  # type: ignore[misc]
    with pytest.raises(AttributeError):
        calendar.anything = 1  # type: ignore[attr-defined]


# ---------------------------------------------------------------- panel.term_arrays


def _probe() -> xt.PanelData:
    frame = pd.DataFrame(
        {
            "id": [2, 1, 1, 2, 1, 3, 3],
            "t": [1, 3, 1, 2, 2, 2, 3],
            "x": [5.0, 3.0, 1.0, 4.0, 2.0, 7.0, 8.0],
            "z": [1.0, np.nan, 3.0, 4.0, 5.0, 6.0, 7.0],
        }
    )
    return xt.xtset(frame, "id", "t")


def test_term_arrays_of_plain_columns_equals_arrays() -> None:
    panel = _probe()
    terms, plain = panel.term_arrays(["x", "z"]), panel.arrays(["x", "z"])
    assert isinstance(terms, PanelArrays)
    np.testing.assert_array_equal(terms.values, plain.values)
    assert terms.columns == plain.columns
    for field in ("offsets", "unit_of_row", "time_index", "source_row", "valid"):
        np.testing.assert_array_equal(getattr(terms, field), getattr(plain, field))
    assert (terms.n_source, terms.delta, terms.contract_version) == (
        plain.n_source,
        plain.delta,
        plain.contract_version,
    )


def test_term_arrays_equals_operate_in_the_order_of_the_arrays() -> None:
    panel = _probe()
    terms = panel.term_arrays(["x", "L.x", "D.x"])
    assert terms.columns == ("x", "L.x", "D.x")
    expected = panel.operate("L.x", "D.x").assign(x=panel.frame["x"])[["x", "L.x", "D.x"]]
    np.testing.assert_array_equal(terms.values, expected.iloc[terms.source_row].to_numpy())


def test_term_arrays_accepts_one_string_and_canonical_names() -> None:
    panel = _probe()
    assert panel.term_arrays("L1.x").columns == ("L.x",)
    assert panel.term_arrays("x").values.shape == (panel.n_obs, 1)


def test_term_arrays_marks_the_rows_a_lag_leaves_missing() -> None:
    arrays = _probe().term_arrays(["L.x"])
    assert arrays.valid.sum() < arrays.n_obs  # each unit loses its first period


# ---------------------------------------------------------------- panel.sorted_column


def _labelled() -> xt.PanelData:
    frame = pd.DataFrame(
        {
            "id": [2, 1, 1, 2, 1, 3, 3],
            "t": [1, 3, 1, 2, 2, 2, 3],
            "x": [5.0, 3.0, 1.0, 4.0, 2.0, 7.0, 8.0],
            "region": ["south", "north", "north", "south", "north", None, None],
            "when": pd.to_datetime(
                ["2001-01-01", "2003-01-01", "2001-01-01", "2002-01-01", "2002-01-01"] + [None] * 2
            ),
            "kind": pd.Categorical(["a", "b", "b", "a", "b", "c", "c"]),
        }
    )
    return xt.xtset(frame, "id", "t", missing_keys="error")


def test_sorted_column_follows_the_rows_of_arrays_for_any_type() -> None:
    panel = _labelled()
    rows = panel.arrays(["x"]).source_row
    for name in ("region", "when", "kind", "x"):
        expected = panel.frame[name].to_numpy()[rows]
        pd.testing.assert_series_equal(
            pd.Series(panel.sorted_column(name)), pd.Series(expected), check_dtype=False
        )
    region = panel.sorted_column("region")
    assert region[:5].tolist() == ["north", "north", "north", "south", "south"]
    assert pd.isna(region[5:]).all()  # missing values are kept as they are


def test_sorted_column_reads_the_key_columns_and_index_levels() -> None:
    panel = _labelled()
    assert panel.sorted_column("id").tolist() == [1, 1, 1, 2, 2, 3, 3]
    indexed = xt.xtset(panel.frame.set_index(["id", "t"]))
    assert indexed.sorted_column("t").tolist() == panel.sorted_column("t").tolist()


def test_sorted_column_is_a_copy_and_an_unknown_name_proposes_close_ones() -> None:
    panel = _labelled()
    taken = panel.sorted_column("x")
    taken[:] = -1.0
    assert (panel.sorted_column("x") != -1.0).all()  # the panel is untouched
    with pytest.raises(DataError) as error:
        panel.sorted_column("regon")
    assert "region" in str(error.value)


# ---------------------------------------------------------------- panel.numeric_names


def test_numeric_names_lists_the_numeric_columns_in_column_order_without_the_keys() -> None:
    frame = pd.DataFrame(
        {
            "region": ["a", "a", "b", "b"],
            "z": [1.0, 2.0, 3.0, 4.0],
            "id": [1, 1, 2, 2],
            "t": [1, 2, 1, 2],
            "flag": [True, False, True, True],
            "n": [1, 2, 3, 4],
            "when": pd.to_datetime(["2001-01-01"] * 4),
            "kind": pd.Categorical(["a", "b", "a", "b"]),
        }
    )
    panel = xt.xtset(frame, "id", "t")
    assert panel.numeric_names == ("z", "flag", "n")  # no keys, no text, date or category
    assert isinstance(panel.numeric_names, tuple)


def test_every_numeric_name_can_be_read_by_arrays_and_agrees_with_variation() -> None:
    panel = _labelled()
    names = panel.numeric_names
    assert names == ("x",)
    assert panel.arrays(list(names)).columns == names
    assert list(panel.variation().table.index) == list(names)


def test_a_panel_with_only_its_keys_has_no_numeric_names() -> None:
    panel = xt.xtset(pd.DataFrame({"id": [1, 1, 2], "t": [1, 2, 1]}), "id", "t")
    assert panel.numeric_names == ()


def test_numeric_names_follow_the_selection_not_a_stale_frame() -> None:
    panel = _probe()
    assert panel.select("x > 0").numeric_names == panel.numeric_names


# ---------------------------------------------------------------- protocols


def test_a_panel_satisfies_extended_panel_like_and_panel_like_is_not_widened() -> None:
    panel = _probe()
    assert isinstance(panel, extend.ExtendedPanelLike)
    assert isinstance(panel, xt.PanelLike)

    class OldStyle:
        id_name, time_name, n_units, n_obs = "id", "t", 1, 1

        def arrays(self, columns: object, *, where: object = None) -> object:
            return None

    assert isinstance(OldStyle(), xt.PanelLike)  # an estimator-side stand-in still qualifies
    assert not isinstance(OldStyle(), extend.ExtendedPanelLike)
