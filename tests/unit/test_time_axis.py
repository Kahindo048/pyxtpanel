import warnings

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from pyxtpanel._config import option_context
from pyxtpanel._errors import PyxtpanelWarning, SpecificationError, TimeAxisError
from pyxtpanel.time.axis import check_span, encode_time
from pyxtpanel.time.formats import freq_from_format, label_period


def ranks(values, **kw):  # type: ignore[no-untyped-def]
    return encode_time(values, **kw).rank.tolist()


# ---------------------------------------------------------------- integers


def test_integer_time_gives_calendar_ranks() -> None:
    enc = encode_time(np.array([1, 2, 3, 5]))
    assert enc.rank.tolist() == [0, 1, 2, 4]
    assert (enc.axis.kind, enc.axis.freq, enc.axis.delta, enc.axis.n_steps) == (
        "integer",
        None,
        1,
        5,
    )


def test_ranks_do_not_depend_on_row_order() -> None:
    assert ranks(np.array([5, 1, 3, 2])) == [4, 0, 2, 1]


def test_float_time_with_integer_values_is_accepted() -> None:
    assert ranks(np.array([2.0, 3.0, 4.0])) == [0, 1, 2]


def test_decimal_time_is_refused_like_stata_451() -> None:
    units = np.array(["a", "a"])
    with pytest.raises(TimeAxisError, match=r"2\.5.*row 1.*\(firm a\)") as info:
        encode_time(np.array([1.0, 2.5]), name="year", units=units, unit_name="firm")
    assert info.value.stata_code == 451 and info.value.key == "time.not_integer"


def test_infinite_time_is_refused() -> None:
    with pytest.raises(TimeAxisError, match="integer values"):
        encode_time(np.array([1.0, np.inf]))


def test_text_time_is_refused() -> None:
    with pytest.raises(TimeAxisError, match="numeric, a date or a period") as info:
        encode_time(np.array(["2020", "2021"]))
    assert info.value.key == "time.not_numeric"


def test_delta_two_ranks_and_off_grid_refusal() -> None:
    assert ranks(np.array([1, 3, 7]), delta=2) == [0, 1, 3]
    with pytest.raises(TimeAxisError, match=r"value 4 \(row 1\).*multiple of 2") as info:
        encode_time(np.array([1, 4, 7]), delta=2)
    assert info.value.stata_code == 451 and info.value.key == "time.not_aligned"


def test_delta_message_in_french() -> None:
    with option_context(lang="fr"), pytest.raises(TimeAxisError, match="n'est pas sur la grille"):
        encode_time(np.array([1, 4]), delta=2)


@pytest.mark.parametrize("bad", [0, -1, 1.5, True, "2"])
def test_invalid_delta(bad) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(SpecificationError) as info:
        encode_time(np.array([1, 2]), delta=bad)
    assert info.value.stata_code == 198 and info.value.key == "time.delta_invalid"


def test_suggests_delta_when_every_value_is_a_multiple() -> None:
    with pytest.warns(PyxtpanelWarning, match="delta=2"):
        encode_time(np.array([2, 4, 6, 8]))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        encode_time(np.array([2, 4, 6, 8]), delta=2)
        encode_time(np.array([2, 3, 4]))


def test_a_single_period_has_one_step() -> None:
    enc = encode_time(np.array([1990, 1990, 1990]))
    assert enc.rank.tolist() == [0, 0, 0] and enc.axis.n_steps == 1


def test_empty_input() -> None:
    enc = encode_time(np.array([], dtype=np.int64))
    assert enc.rank.size == 0 and enc.axis.n_steps == 0


def test_freq_with_plain_integers_is_refused() -> None:
    with pytest.raises(SpecificationError, match="plain integers") as info:
        encode_time(np.array([1, 2]), freq="M")
    assert info.value.key == "time.freq_without_dates"


def test_unknown_freq() -> None:
    with pytest.raises(SpecificationError, match="Unknown frequency"):
        encode_time(pd.to_datetime(["2020-01-01"]), freq="X")


# ---------------------------------------------------------------- dates


def dates(*text: str) -> pd.Series:
    return pd.Series(pd.to_datetime(list(text)))


@pytest.mark.parametrize(
    ("text", "freq", "expected"),
    [
        (("2019-12-01", "2020-01-01", "2020-03-01"), "M", [0, 1, 3]),
        (("2019-10-01", "2020-01-01", "2020-07-01"), "Q", [0, 1, 3]),
        (("2018-01-01", "2018-07-01", "2019-07-01"), "H", [0, 1, 3]),
        (("2018-01-01", "2019-01-01", "2021-01-01"), "Y", [0, 1, 3]),
        (("2020-01-06", "2020-01-13", "2020-01-27"), "W", [0, 1, 3]),
        (("2020-02-27", "2020-02-28", "2020-03-01"), "D", [0, 1, 3]),
    ],
)
def test_frequency_is_inferred_and_ranks_follow_the_calendar(text, freq, expected) -> None:  # type: ignore[no-untyped-def]
    enc = encode_time(dates(*text))
    assert enc.axis.freq == freq and enc.rank.tolist() == expected
    assert enc.axis.kind == "datetime"


def test_month_end_dates_are_monthly() -> None:
    enc = encode_time(dates("2020-01-31", "2020-02-29", "2020-04-30"))
    assert enc.axis.freq == "M" and enc.rank.tolist() == [0, 1, 3]


def test_explicit_freq_overrides_inference_and_accepts_any_day() -> None:
    enc = encode_time(dates("2020-01-15", "2020-02-03", "2020-04-30"), freq="M")
    assert enc.rank.tolist() == [0, 1, 3]


def test_quarter_with_delta() -> None:
    enc = encode_time(dates("2020-01-01", "2020-07-01", "2021-07-01"), freq="Q", delta=2)
    assert (
        enc.axis.freq == "Q"
        and enc.rank.tolist() == [0, 1, 3]
        and enc.axis.delta_text == "2 quarters"
    )
    with pytest.raises(TimeAxisError, match="not on the calendar grid"):
        encode_time(dates("2020-01-01", "2020-04-01"), freq="Q", delta=2)


def test_period_dtype_carries_its_frequency() -> None:
    periods = pd.Series(pd.period_range("2020Q1", periods=4, freq="Q"))
    enc = encode_time(periods)
    assert enc.axis.kind == "period" and enc.axis.freq == "Q" and enc.rank.tolist() == [0, 1, 2, 3]
    assert list(enc.axis.values(enc.rank)) == list(periods)


def test_intraday_times_are_refused_for_daily_calendar() -> None:
    with pytest.raises(TimeAxisError, match="time of day"):
        encode_time(dates("2020-01-01 10:00", "2020-01-02 10:00"), freq="D")


def test_missing_date_is_refused() -> None:
    series = pd.Series(pd.to_datetime(["2020-01-01", None]))
    with pytest.raises(TimeAxisError, match="missing"):
        encode_time(series)


def test_business_days_are_flagged() -> None:
    days = pd.bdate_range("2020-01-01", periods=60)
    with pytest.warns(PyxtpanelWarning, match="Saturday or a Sunday"):
        enc = encode_time(pd.Series(days))
    assert enc.axis.freq == "D" and enc.rank.max() > 59  # weekends are gaps
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        encode_time(pd.Series(pd.date_range("2020-01-01", periods=60)))


def test_values_and_labels_roundtrip() -> None:
    series = dates("2019-12-01", "2020-01-01", "2020-04-01")
    enc = encode_time(series)
    assert enc.axis.labels(enc.rank) == ["2019m12", "2020m1", "2020m4"]
    assert list(enc.axis.values(enc.rank)) == list(pd.DatetimeIndex(series))
    quarterly = encode_time(dates("2019-10-01", "2020-01-01"))
    assert quarterly.axis.labels(quarterly.rank) == ["2019q4", "2020q1"]
    daily = encode_time(dates("2020-02-27", "2020-03-01"))
    assert daily.axis.labels(daily.rank) == ["27feb2020", "01mar2020"]
    yearly = encode_time(dates("2019-01-01", "2021-01-01"))
    assert yearly.axis.labels(yearly.rank) == ["2019", "2021"]


# ---------------------------------------------------------------- Stata-coded integers


def test_stata_coded_months_and_quarters() -> None:
    enc = encode_time(np.array([720, 721, 723]), format="%tm")  # 2020m1, 2020m2, 2020m4
    assert enc.axis.kind == "stata" and enc.axis.freq == "M" and enc.rank.tolist() == [0, 1, 3]
    assert enc.axis.labels(enc.rank) == ["2020m1", "2020m2", "2020m4"]
    assert enc.axis.values(enc.rank).tolist() == [720, 721, 723]
    quarters = encode_time(np.array([240, 243]), format="%tq")
    assert quarters.axis.labels(quarters.rank) == ["2020q1", "2020q4"]
    years = encode_time(np.array([1987, 1990]), format="%ty")
    assert years.rank.tolist() == [0, 3] and years.axis.values(years.rank).tolist() == [1987, 1990]
    days = encode_time(np.array([0, 1, 22000]), format="%td")  # 01jan1960 is day 0
    assert days.axis.labels(days.rank)[0] == "01jan1960"


@pytest.mark.parametrize("fmt", ["%tw", "%tc", "%tC", "%9.0g", "td"])
def test_unsupported_stata_formats(fmt) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(SpecificationError, match="not supported") as info:
        encode_time(np.array([1, 2]), format=fmt)
    assert info.value.stata_code == 198


def test_format_parsing_accepts_display_variants() -> None:
    assert freq_from_format("%tqCCYY!Qq") == "Q" and freq_from_format("%-tm") == "M"
    assert label_period(0, "Y") == "1970"


# ---------------------------------------------------------------- span guard


def test_check_span_refuses_overflowing_keys() -> None:
    axis = encode_time(np.array([0, 10**12])).axis
    check_span(axis, 1000)
    with pytest.raises(TimeAxisError, match="too large"):
        check_span(axis, 10**8)


# ---------------------------------------------------------------- properties


@given(
    start=st.integers(-1000, 1000),
    delta=st.integers(1, 7),
    steps=st.lists(st.integers(0, 40), min_size=1, max_size=30),
)
def test_rank_times_delta_recovers_the_value(start: int, delta: int, steps: list[int]) -> None:
    values = np.array([start + delta * s for s in steps], dtype=np.int64)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PyxtpanelWarning)
        enc = encode_time(values, delta=delta)
    assert (enc.axis.values(enc.rank).to_numpy() == values).all()
    assert (np.diff(np.sort(enc.rank)) >= 0).all()
    assert enc.rank.min() == 0 and enc.axis.n_steps == enc.rank.max() + 1


@given(months=st.lists(st.integers(0, 600), min_size=1, max_size=30))
def test_monthly_dates_rank_equals_month_difference(months: list[int]) -> None:
    stamps = pd.Series(pd.to_datetime([f"{1970 + m // 12}-{m % 12 + 1:02d}-01" for m in months]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PyxtpanelWarning)
        enc = encode_time(stamps, freq="M")
    assert enc.rank.tolist() == [m - min(months) for m in months]
