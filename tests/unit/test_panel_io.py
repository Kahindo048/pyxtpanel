from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.panel.base import PanelBase
from pyxtpanel.panel.io import (
    _hex_delta,
    _stata_delta,
    inject_xtset,
    read_stata,
    read_table,
    stata_xtset,
)

INVEST = Path(__file__).resolve().parents[2] / "validation" / "data" / "invest.dta"


def make() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": [1, 1, 1, 2, 2],
            "t": [2000, 2001, 2003, 2000, 2001],
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
        }
    )


# ---------------------------------------------------------------- hexadecimal delta


@pytest.mark.parametrize("delta", [1, 2, 3, 4, 7, 12, 100, 1000])
def test_hex_delta_roundtrips(delta: int) -> None:
    assert _stata_delta(_hex_delta(delta)) == delta


def test_hex_delta_matches_the_text_stata_writes() -> None:
    assert _hex_delta(1) == "+1.0000000000000X+000"
    assert _hex_delta(2) == "+1.0000000000000X+001"
    assert _hex_delta(3) == "+1.8000000000000X+001"


def test_unparsable_delta() -> None:
    assert _stata_delta("garbage") is None


# ---------------------------------------------------------------- reading xtset from a .dta


@pytest.mark.skipif(not INVEST.exists(), reason="validation data not present")
def test_reads_stata_made_file_without_xtset() -> None:
    assert stata_xtset(INVEST) == stata_xtset(INVEST).__class__(None, None, None)
    frame = read_stata(INVEST)
    assert frame.attrs["id"] is None and frame.attrs["time"] is None
    with pytest.raises(SpecificationError, match="no xtset declaration") as info:
        PanelBase.from_stata(INVEST)
    assert info.value.stata_code == 459
    assert PanelBase.from_stata(INVEST, "firm", "year").n_units == 565


def test_unreadable_or_foreign_files_give_empty_info(tmp_path: Path) -> None:
    junk = tmp_path / "junk.dta"
    junk.write_bytes(b"not a stata file")
    assert stata_xtset(junk).panel is None
    assert stata_xtset(tmp_path / "missing.dta").panel is None


# ---------------------------------------------------------------- writing and reading back


def test_write_then_read_keeps_the_declaration(tmp_path: Path) -> None:
    path = tmp_path / "panel.dta"
    panel = PanelBase.from_frame(make(), "id", "t")
    panel.to_stata(path)
    info = stata_xtset(path)
    assert (info.panel, info.time, info.delta) == ("id", "t", 1.0)
    again = PanelBase.from_stata(path)
    assert (again.id_name, again.time_name, again.n_obs) == ("id", "t", 5)
    assert again.has_gaps and again.balanced == "unbalanced"
    pd.testing.assert_frame_equal(
        again.to_frame(order="sorted").reset_index(drop=True),
        panel.to_frame(order="sorted").reset_index(drop=True),
        check_dtype=False,
    )


def test_delta_is_recorded_and_read_back(tmp_path: Path) -> None:
    path = tmp_path / "delta.dta"
    frame = pd.DataFrame({"id": [1, 1, 1], "t": [1, 3, 7]})
    PanelBase.from_frame(frame, "id", "t", delta=2).to_stata(path)
    assert stata_xtset(path).delta == 2.0
    again = PanelBase.from_stata(path)
    assert again.delta == 2 and again.n_gaps == 1


def test_dates_are_written_with_their_stata_format(tmp_path: Path) -> None:
    path = tmp_path / "dates.dta"
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 2],
            "q": pd.to_datetime(["2020-01-01", "2020-04-01", "2020-10-01", "2020-01-01"]),
            "x": [1.0, 2.0, 3.0, 4.0],
        }
    )
    PanelBase.from_frame(frame, "id", "q").to_stata(path)
    assert read_stata(path).attrs["time_format"] == "%tq"
    again = PanelBase.from_stata(path)
    assert again.axis.freq == "Q" and again.n_gaps == 1 and again.n_obs == 4
    assert again.t_min == pd.Timestamp("2020-01-01")


def test_writing_refuses_text_ids_and_weekly_data(tmp_path: Path) -> None:
    text = pd.DataFrame({"id": ["a", "a"], "t": [1, 2]})
    with pytest.raises(SpecificationError, match="numeric panel variable"):
        PanelBase.from_frame(text, "id", "t").to_stata(tmp_path / "a.dta")
    weekly = pd.DataFrame(
        {"id": [1, 1, 1], "d": pd.to_datetime(["2020-01-06", "2020-01-13", "2020-01-20"])}
    )
    with pytest.raises(SpecificationError, match="Weekly"):
        PanelBase.from_frame(weekly, "id", "d").to_stata(tmp_path / "w.dta")


def test_stata_weeks_survive_a_round_trip_through_a_dta_file(tmp_path: Path) -> None:
    """``%tw`` is written as Stata's weekly format and read back as the same panel."""
    frame = pd.DataFrame(
        {
            "id": [1, 1, 1, 2, 2],
            "w": [2078, 2079, 2081, 2079, 2080],  # 1999w51 ... 2000w2, across the new year
            "x": [1.0, 2.0, 4.0, 5.0, 6.0],
        }
    )
    panel = PanelBase.from_frame(frame, "id", "w", format="%tw")
    path = tmp_path / "weeks.dta"
    panel.to_stata(path)
    assert read_stata(path).attrs["time_format"] == "%tw"
    again = PanelBase.from_stata(path)
    assert again.axis.freq == "SW" and again.n_gaps == panel.n_gaps == 1
    assert again.axis.labels(again.index.rank)[:3] == ["1999w51", "1999w52", "2000w2"]
    np.testing.assert_array_equal(again.index.rank, panel.index.rank)


def test_inject_xtset_into_a_foreign_file_is_refused(tmp_path: Path) -> None:
    bad = tmp_path / "bad.dta"
    bad.write_bytes(b"nope")
    with pytest.raises(DataError):
        inject_xtset(bad, "id", "t", 1)


def test_data_survive_the_injection(tmp_path: Path) -> None:
    path = tmp_path / "plain.dta"
    frame = make()
    frame.to_stata(path, write_index=False, version=118)
    before = pd.read_stata(path)
    inject_xtset(path, "id", "t", 3)
    pd.testing.assert_frame_equal(pd.read_stata(path), before)
    assert stata_xtset(path).delta == 3.0


def test_override_ids_and_non_integer_recorded_delta(tmp_path: Path) -> None:
    path = tmp_path / "o.dta"
    PanelBase.from_frame(make(), "id", "t").to_stata(path)
    other = PanelBase.from_stata(path, id="t", time="id", missing_keys="drop")
    assert other.id_name == "t"
    plain = tmp_path / "half.dta"
    make().to_stata(plain, write_index=False, version=118)
    inject_xtset(plain, "id", "t", 1)
    raw = plain.read_bytes().replace(_hex_delta(1).encode(), b"+1.0000000000000X-001")
    plain.write_bytes(raw)
    with pytest.raises(DataError, match="not a whole number"):
        PanelBase.from_stata(plain)
    assert PanelBase.from_stata(plain, delta=1).n_obs == 5


# ---------------------------------------------------------------- other formats


def test_read_table_csv_and_unknown(tmp_path: Path) -> None:
    csv = tmp_path / "a.csv"
    make().to_csv(csv, index=False)
    assert read_table(csv).shape == (5, 3)
    with pytest.raises(DataError, match="Unsupported file type"):
        read_table(tmp_path / "a.xlsx")
    stata = tmp_path / "a.dta"
    make().to_stata(stata, write_index=False)
    assert read_table(stata).shape == (5, 3)


def test_read_parquet_when_available(tmp_path: Path) -> None:
    pytest.importorskip("pyarrow")
    path = tmp_path / "a.parquet"
    make().to_parquet(path)
    assert np.array_equal(read_table(path)["x"], make()["x"])
