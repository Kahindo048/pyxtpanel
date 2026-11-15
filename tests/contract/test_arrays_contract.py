"""Contract test of the array surface used by the estimators (D.3.3).

Estimation packages depend on this layout; the test fails on any change that would break them.
A deliberate change must raise ``CONTRACT_VERSION`` and update the snapshot below.
"""

import dataclasses
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

from pyxtpanel import PyxtpanelWarning, xtset
from pyxtpanel._errors import DataError
from pyxtpanel.contract import CONTRACT_VERSION, PanelArrays
from pyxtpanel.ops.diff import fd_equations

FIELDS_V1 = {
    "values",
    "columns",
    "offsets",
    "unit_labels",
    "unit_of_row",
    "time_index",
    "source_row",
    "valid",
    "n_source",
    "delta",
    "contract_version",
}
METHODS_V1 = {"unit_slices", "lag_index", "reindex_back", "as_cube"}


def make():  # type: ignore[no-untyped-def]
    frame = pd.DataFrame(
        {
            "id": [2, 1, 1, 2, 1, 3],
            "t": [3, 1, 5, 1, 3, 5],
            "x": [1.0, 2, 3, 4, 5, 6],
            "z": [np.nan, 1, 2, 3, 4, 5],
        },
        index=list("abcdef"),
    )
    return frame, xtset(frame, "id", "t", delta=2)


def test_the_version_and_the_field_snapshot_are_pinned() -> None:
    assert CONTRACT_VERSION == 1
    assert {f.name for f in dataclasses.fields(PanelArrays)} == FIELDS_V1
    assert set(dir(PanelArrays)) >= METHODS_V1
    assert all(f.kw_only for f in dataclasses.fields(PanelArrays))  # fields can be added safely


def test_layout_dtypes_and_order() -> None:
    frame, panel = make()
    arrays = panel.arrays(["x", "z"])
    assert arrays.contract_version == CONTRACT_VERSION
    assert arrays.values.dtype == np.float64 and arrays.values.flags.c_contiguous
    assert arrays.values.shape == (6, 2) and arrays.columns == ("x", "z")
    for name in ("offsets", "unit_of_row", "time_index", "source_row"):
        array = getattr(arrays, name)
        assert array.dtype == np.int64 and not array.flags.writeable, name
    assert arrays.valid.dtype == np.bool_
    assert arrays.offsets.tolist() == [0, 3, 5, 6] and arrays.unit_labels.tolist() == [1, 2, 3]
    # sorted by unit then calendar position, strictly increasing within a unit
    assert arrays.time_index.tolist() == [0, 1, 2, 0, 1, 2]
    for part in arrays.unit_slices():
        assert np.all(np.diff(arrays.time_index[part]) > 0)
    assert np.all(np.diff(arrays.unit_of_row) >= 0)
    assert arrays.valid.tolist() == [True, True, True, True, False, True]
    assert arrays.n_obs == 6 and arrays.n_units == 3 and arrays.counts.tolist() == [3, 2, 1]


def test_time_index_is_the_calendar_position_whatever_delta() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1], "t": [10, 14, 18], "x": [1.0, 2, 3]})
    arrays = xtset(frame, "id", "t", delta=4).arrays(["x"])
    assert arrays.time_index.tolist() == [0, 1, 2] and arrays.delta == 4
    dates = pd.DataFrame(
        {"id": [1, 1], "t": pd.to_datetime(["2020-01-01", "2021-01-01"]), "x": [1.0, 2]}
    )
    assert xtset(dates, "id", "t", freq="Y").arrays(["x"]).time_index.tolist() == [0, 1]


def test_index_arrays_are_shared_and_values_are_a_private_copy() -> None:
    frame, panel = make()
    arrays = panel.arrays(["x"])
    assert np.shares_memory(arrays.offsets, panel.index.offsets)
    assert np.shares_memory(arrays.time_index, panel.index.rank)
    arrays.values[0, 0] = -99.0  # the caller owns the values
    assert panel.arrays(["x"]).values[0, 0] != -99.0 and frame["x"].iloc[0] != -99.0
    with pytest.raises(ValueError, match="read-only"):
        arrays.offsets[0] = 5


def test_reindex_back_restores_the_original_order_and_pads_missing_rows() -> None:
    frame, panel = make()
    arrays = panel.arrays(["x"])
    back = arrays.reindex_back(arrays.values[:, 0])
    assert back.tolist() == frame["x"].tolist()
    sub = panel.arrays(["x"], where="id != 1")
    out = sub.reindex_back(sub.values[:, 0])
    assert out.shape == (6,) and np.isnan(out).sum() == 3 and out[0] == 1.0
    two = arrays.reindex_back(np.column_stack([arrays.values[:, 0]] * 2))
    assert two.shape == (6, 2)
    with pytest.raises(DataError, match="expected 6 rows"):
        arrays.reindex_back(np.zeros(5))


def test_lag_index_follows_the_calendar_not_the_previous_row() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1, 2, 2], "t": [1, 2, 4, 1, 2], "x": [1.0, 2, 3, 4, 5]})
    arrays = xtset(frame, "id", "t").arrays(["x"])
    assert arrays.lag_index(1).tolist() == [-1, 0, -1, -1, 3]
    assert arrays.lag_index(2).tolist() == [-1, -1, 1, -1, -1]
    assert arrays.lag_index(-1).tolist() == [1, -1, -1, 4, -1]
    assert arrays.lag_index(0).tolist() == [0, 1, 2, 3, 4]


def test_unit_slices_cover_the_rows_exactly_once() -> None:
    _, panel = make()
    arrays = panel.arrays(["x"])
    rows = np.concatenate([np.arange(6)[s] for s in arrays.unit_slices()])
    assert rows.tolist() == list(range(6))


def test_as_cube_only_for_strongly_balanced_panels() -> None:
    frame = pd.DataFrame(
        {"id": [1, 1, 2, 2], "t": [1, 3, 1, 3], "x": [1.0, 2, 3, 4], "z": [5.0, 6, 7, 8]}
    )
    cube = xtset(frame, "id", "t").arrays(["x", "z"]).as_cube()
    assert cube.shape == (2, 2, 2) and cube[1, 1].tolist() == [4.0, 8.0]
    _, unbalanced = make()
    with pytest.raises(DataError, match="strongly balanced"):
        unbalanced.arrays(["x"]).as_cube()


def test_infinite_values_become_missing_with_a_warning() -> None:
    frame = pd.DataFrame({"id": [1, 1], "t": [1, 2], "x": [1.0, np.inf]})
    panel = xtset(frame, "id", "t")
    with pytest.warns(PyxtpanelWarning, match="1 infinite"):
        arrays = panel.arrays(["x"])
    assert np.isnan(arrays.values[1, 0]) and not arrays.valid[1]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        panel.select(complete=["id"]).arrays(["t"])  # no infinity, no warning


def test_first_difference_equations_pair_each_row_with_its_previous_period() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1, 2, 2], "t": [1, 2, 4, 1, 2], "x": [1.0, 2, 3, 4, 5]})
    panel = xtset(frame, "id", "t")
    rows, previous = fd_equations(panel.index)
    assert rows.tolist() == [1, 4] and previous.tolist() == [0, 3]


def test_importing_the_contract_loads_no_compiled_kernel() -> None:
    code = (
        "import sys; import pyxtpanel.contract; "
        "sys.exit(int(any(m == 'numba' or m.startswith('numba.') for m in sys.modules)))"
    )
    assert subprocess.run([sys.executable, "-c", code], check=False).returncode == 0
