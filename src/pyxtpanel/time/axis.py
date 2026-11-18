"""The time axis of a panel: every period becomes an integer rank on a regular calendar.

A panel is observed on a calendar, not on a list of rows.  :func:`encode_time` turns a time
column into integer ranks ``0, 1, 2, ...`` such that two periods that are ``k`` steps apart on
the calendar have ranks that differ by ``k``.  A lag is then a subtraction on ranks
(``t - 1``), whatever the unit of the time variable: an integer with a step ``delta``, dates
with a frequency (day, week, Stata week, month, quarter, half-year, year), ``pandas.Period``
values, or Stata-coded integers with a ``%t`` format.

The axis is *global*: every unit shares the same calendar.  This is stricter than Stata, which
only requires, within each unit, differences that are at least ``delta``.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from functools import reduce
from math import gcd
from typing import Any, Literal, cast

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import PyxtpanelWarning, SpecificationError, TimeAxisError
from pyxtpanel._types import IntArray
from pyxtpanel.time.formats import STATA_OFFSET, Freq, freq_from_format, label_period

FREQS: tuple[Freq, ...] = ("D", "W", "SW", "M", "Q", "H", "Y")
Kind = Literal["integer", "datetime", "period", "stata"]

#: Largest ``n_units * n_steps`` for which the lookup keys of the calendar fit in an int64.
MAX_KEY = 2**62

_FREQ_NAMES = {
    "D": "day",
    "W": "week",
    "SW": "week",
    "M": "month",
    "Q": "quarter",
    "H": "half-year",
    "Y": "year",
}

_MAX_EXAMPLES = 3

_i18n.register(
    "time",
    {
        "time.not_numeric": {
            "en": "The time variable '{name}' must be numeric, a date or a period (got {dtype}).",
            "fr": "La variable de temps « {name} » doit être numérique, une date ou une période "
            "(reçu {dtype}).",
        },
        "time.not_integer": {
            "en": "The time variable '{name}' must contain only integer values "
            "(found {example} in row {row}{unit}).",
            "fr": "La variable de temps « {name} » ne doit contenir que des entiers "
            "(trouvé {example} à la ligne {row}{unit}).",
        },
        "time.missing": {
            "en": "The time variable '{name}' has {n} missing value(s) (first rows: {rows}).",
            "fr": "La variable de temps « {name} » a {n} valeur(s) manquante(s) "
            "(premières lignes : {rows}).",
        },
        "time.delta_invalid": {
            "en": "delta must be a positive integer, got {delta!r}.",
            "fr": "delta doit être un entier strictement positif, reçu {delta!r}.",
        },
        "time.freq_invalid": {
            "en": "Unknown frequency {freq!r}; use one of {allowed}.",
            "fr": "Fréquence inconnue {freq!r} ; utilisez l'une de {allowed}.",
        },
        "time.freq_without_dates": {
            "en": "freq={freq!r} needs dates or periods, but '{name}' holds plain integers "
            "(for Stata-coded values, pass format='%t..').",
            "fr": "freq={freq!r} suppose des dates ou des périodes, mais « {name} » contient des "
            "entiers (pour des valeurs codées à la Stata, passez format='%t..').",
        },
        "time.intraday": {
            "en": "Dates with a time of day cannot be used with a daily frequency "
            "(found {example}); normalise them to midnight or give a coarser freq.",
            "fr": "Des dates avec une heure ne peuvent pas servir à une fréquence quotidienne "
            "(trouvé {example}) ; ramenez-les à minuit ou donnez une fréquence plus grossière.",
        },
        "time.not_aligned": {
            "en": "Time value {example} (row {row}{unit}) is not on the calendar grid: with "
            "delta={delta} every value must be a multiple of {delta} steps after the first "
            "period {first}. Check delta= or the data.",
            "fr": "La valeur de temps {example} (ligne {row}{unit}) n'est pas sur la grille du "
            "calendrier : avec delta={delta}, toute valeur doit être à un multiple de {delta} "
            "pas de la première période {first}. Vérifiez delta= ou les données.",
        },
        "time.suggest_delta": {
            "en": "All time values are {delta} steps apart or a multiple of it; with delta=1 "
            "every period looks like a gap. Pass delta={delta} if the panel is observed every "
            "{delta} periods.",
            "fr": "Toutes les valeurs de temps sont espacées de {delta} pas ou d'un multiple ; "
            "avec delta=1 chaque période ressemble à un trou. Passez delta={delta} si le panel "
            "est observé toutes les {delta} périodes.",
        },
        "time.business_days": {
            "en": "No date falls on a Saturday or a Sunday: with a daily calendar every weekend "
            "is a gap, so L.x is missing on Mondays. Use an integer day counter if you want "
            "trading days to be consecutive.",
            "fr": "Aucune date ne tombe un samedi ou un dimanche : avec un calendrier quotidien "
            "chaque week-end est un trou, donc L.x est manquant les lundis. Utilisez un compteur "
            "entier de jours si les jours ouvrés doivent être consécutifs.",
        },
        "time.span_too_large": {
            "en": "The time axis spans {n_steps} steps for {n_units} units: too large to index. "
            "Use a smaller unit of time or a larger delta.",
            "fr": "L'axe du temps couvre {n_steps} pas pour {n_units} unités : trop grand pour "
            "être indexé. Utilisez une unité de temps plus grande ou un delta plus grand.",
        },
        "time.format_unsupported": {
            "en": "Stata format {fmt!r} is not supported (supported: %td %tw %tm %tq %th %ty).",
            "fr": "Le format Stata {fmt!r} n'est pas géré (gérés : %td %tw %tm %tq %th %ty).",
        },
    },
)


@dataclass(frozen=True)
class TimeAxis:
    """A regular calendar shared by all the units of a panel.

    Attributes
    ----------
    kind
        ``"integer"`` (plain integers), ``"datetime"`` (``datetime64``), ``"period"``
        (``pandas.Period``) or ``"stata"`` (integers coded as Stata does for a ``%t`` format).
    freq
        ``None`` for plain integers, otherwise the calendar unit.
    delta
        Number of calendar units between two consecutive steps (Stata's ``delta()``).
    ordinal_min
        Calendar ordinal of the first period (rank 0).  An ordinal counts calendar units from
        1970 (days, months, quarters, half-years, years; weeks from the Monday of that week;
        Stata weeks, 52 a year, from 1970w1),
        or is the integer itself for ``"integer"``.
    n_steps
        Number of steps from the first to the last period, inclusive (``max rank + 1``).
    stata_format
        The ``%t`` format when ``kind`` is ``"stata"``.
    """

    kind: Kind
    freq: Freq | None
    delta: int
    ordinal_min: int
    n_steps: int
    stata_format: str | None = None

    @property
    def span(self) -> int:
        """Number of periods between the first and the last, inclusive, in steps of delta."""
        return self.n_steps

    def ordinal(self, rank: IntArray | int) -> IntArray:
        """Calendar ordinal of the periods at ``rank``."""
        return np.atleast_1d(self.ordinal_min + np.asarray(rank, dtype=np.int64) * self.delta)

    def values(self, rank: IntArray | int) -> pd.Index:
        """Time values (in the type of the input data) of the periods at ``rank``.

        Dates are the first day of the period.
        """
        ordinal = self.ordinal(rank)
        if self.kind == "integer":
            return pd.Index(ordinal.astype(np.int64))
        assert self.freq is not None
        if self.kind == "stata":
            return pd.Index((ordinal - STATA_OFFSET[self.freq]).astype(np.int64))
        starts = pd.DatetimeIndex(_start_dates(np.atleast_1d(ordinal), self.freq))
        if self.kind == "period":
            return cast("pd.Index", pd.PeriodIndex(starts, freq=_PANDAS_FREQ[self.freq]))
        return cast("pd.Index", starts)

    def stata_values(self, rank: IntArray | int) -> IntArray:
        """Value Stata stores for the periods at ``rank`` (its internal integer for the format).

        For dates this is the count Stata keeps (for example 240 for 2020q1 in ``%tq``); for
        plain integers it is the integer itself.
        """
        ordinal = self.ordinal(rank)
        if self.freq is None:
            return ordinal
        return np.asarray(ordinal - STATA_OFFSET[self.freq], dtype=np.int64)

    def datetimes(self, rank: IntArray | int) -> pd.DatetimeIndex:
        """First day of the periods at ``rank`` as timestamps (dated axes only)."""
        assert self.freq is not None
        return pd.DatetimeIndex(_start_dates(np.atleast_1d(self.ordinal(rank)), self.freq))

    def labels(self, rank: IntArray | int) -> list[str]:
        """Printable labels (``2020q1``, ``1987``, ``01jan2020``) of the periods at ``rank``."""
        ordinal = np.atleast_1d(self.ordinal(rank))
        if self.freq is None:
            return [str(int(v)) for v in ordinal]
        return [label_period(int(v), self.freq) for v in ordinal]

    def rank_of(self, values: pd.Series | np.ndarray | pd.Index | Sequence[Any]) -> IntArray:
        """Rank on this axis of time values given in the type of the input data.

        Ranks may be negative or beyond ``n_steps - 1`` for periods outside the observed range.

        Raises
        ------
        TimeAxisError
            If a value is missing, not numeric/date, or not on the calendar grid.
        """
        arr = np.asarray(values)
        if self.kind in ("datetime", "period"):
            days, _, _ = _dates_from(values)
            if np.isnat(days).any():
                raise TimeAxisError(
                    _i18n.message(
                        "time.missing", name="time", n=int(np.isnat(days).sum()), rows=[]
                    ),
                    key="time.missing",
                )
            assert self.freq is not None
            ordinal = _ordinals_from_dates(days, self.freq)
        elif self.kind == "stata":
            assert self.freq is not None
            ordinal = _integers(arr, "time", None, "id", None) + STATA_OFFSET[self.freq]
        else:
            ordinal = _integers(arr, "time", None, "id", None)
        offset = ordinal - self.ordinal_min
        off_grid = offset % self.delta != 0
        if off_grid.any():
            row = int(np.flatnonzero(off_grid)[0])
            raise TimeAxisError(
                _i18n.message(
                    "time.not_aligned",
                    example=_example(values, row),
                    row=row,
                    unit="",
                    delta=self.delta,
                    first=self.labels(0)[0],
                ),
                stata_code=451,
                key="time.not_aligned",
            )
        return np.asarray(offset // self.delta, dtype=np.int64)

    @property
    def delta_text(self) -> str:
        """Step in words, as Stata prints it (``1 unit``, ``2 units``, ``1 quarter``)."""
        unit = "unit" if self.freq is None else _FREQ_NAMES[self.freq]
        return f"{self.delta} {unit}{'s' if self.delta > 1 else ''}"


_PANDAS_FREQ = {"D": "D", "W": "W-SUN", "M": "M", "Q": "Q-DEC", "Y": "Y-DEC"}


@dataclass(frozen=True)
class EncodedTime:
    """Result of :func:`encode_time`: the axis and the rank of each row."""

    axis: TimeAxis
    rank: IntArray


def _check_delta(delta: int | None) -> int:
    if delta is None:
        return 1
    if isinstance(delta, bool) or not isinstance(delta, (int, np.integer)) or delta < 1:
        raise SpecificationError(
            _i18n.message("time.delta_invalid", delta=delta),
            stata_code=198,
            key="time.delta_invalid",
        )
    return int(delta)


def _check_freq(freq: str | None) -> Freq | None:
    if freq is None:
        return None
    if freq not in FREQS:
        raise SpecificationError(
            _i18n.message("time.freq_invalid", freq=freq, allowed=list(FREQS)),
            stata_code=198,
            key="time.freq_invalid",
        )
    return freq


def _row(position: int, positions: IntArray | None) -> int:
    """Row of the user's frame that holds the ``position``-th encoded value."""
    return position if positions is None else int(positions[position])


def _where(row: int, units: np.ndarray | None, unit_name: str) -> str:
    if units is None or row >= len(units):
        return ""
    return f" ({unit_name} {units[row]})"


def _ordinals_from_dates(days: np.ndarray, freq: Freq) -> IntArray:
    """Calendar ordinals of ``datetime64[D]`` values (the period that contains each date)."""
    if freq == "D":
        return days.astype(np.int64)
    if freq == "W":
        return (days.astype(np.int64) + 3) // 7  # weeks start on Monday; 1970-01-01 is Thursday
    if freq == "SW":  # Stata's week: 52 a year, the 52nd takes the last 8 or 9 days
        years = days.astype("datetime64[Y]")
        day_of_year = (days - years.astype("datetime64[D]")).astype(np.int64)
        return years.astype(np.int64) * 52 + np.minimum(day_of_year // 7, 51)
    months = days.astype("datetime64[M]").astype(np.int64)
    return months // {"M": 1, "Q": 3, "H": 6, "Y": 12}[freq]


def _start_dates(ordinal: IntArray, freq: Freq) -> np.ndarray:
    """First day of the period of each calendar ordinal, as ``datetime64[ns]``."""
    if freq == "D":
        return ordinal.astype("datetime64[D]").astype("datetime64[ns]")
    if freq == "W":
        return ((ordinal * 7 - 3).astype("datetime64[D]")).astype("datetime64[ns]")
    if freq == "SW":
        years = (ordinal // 52).astype("datetime64[Y]").astype("datetime64[D]")
        return (years + (ordinal % 52 * 7).astype("timedelta64[D]")).astype("datetime64[ns]")
    months = ordinal * {"M": 1, "Q": 3, "H": 6, "Y": 12}[freq]
    return months.astype("datetime64[M]").astype("datetime64[ns]")


def _aligned(days: np.ndarray, freq: Freq) -> bool:
    """Whether all dates are period starts (or all period ends), for inferring a frequency."""
    months = days.astype("datetime64[M]")
    first = months.astype("datetime64[D]")
    last = (months + np.timedelta64(1, "M")).astype("datetime64[D]") - np.timedelta64(1, "D")
    month_number = (months.astype(np.int64) % 12) + 1
    if freq == "W":
        unique = np.unique(days)
        return unique.size > 1 and bool(np.all(np.mod(np.diff(unique.astype(np.int64)), 7) == 0))
    allowed = {
        "M": np.ones_like(month_number, dtype=bool),
        "Q": month_number % 3 == 1,
        "H": month_number % 6 == 1,
        "Y": month_number == 1,
    }
    allowed_end = {
        "M": np.ones_like(month_number, dtype=bool),
        "Q": month_number % 3 == 0,
        "H": month_number % 6 == 0,
        "Y": month_number == 12,
    }
    starts = bool(np.all((days == first) & allowed[freq]))
    ends = bool(np.all((days == last) & allowed_end[freq]))
    return starts or ends


def _infer_freq(days: np.ndarray) -> Freq:
    for freq in ("Y", "H", "Q", "M", "W"):
        if _aligned(days, freq):
            return freq
    return "D"


def _dates_from(
    values: pd.Series | np.ndarray | pd.Index | Sequence[Any],
) -> tuple[np.ndarray, Kind, Freq | None]:
    """Return ``datetime64[D]`` days, the kind and the frequency carried by the type."""
    series = pd.Series(values)
    if isinstance(series.dtype, pd.PeriodDtype):
        carried = _freq_of_period(_period_freqstr(series.dtype))
        return series.dt.start_time.to_numpy().astype("datetime64[D]"), "period", carried
    raw = series.to_numpy(dtype="datetime64[ns]")
    days = raw.astype("datetime64[D]")
    return days, "datetime", None


def _period_freqstr(dtype: Any) -> str:
    return str(dtype.freq.freqstr)


def _freq_of_period(freqstr: str) -> Freq:
    head = freqstr.split("-")[0].upper().lstrip("0123456789")
    mapping = {
        "D": "D",
        "W": "W",
        "M": "M",
        "Q": "Q",
        "A": "Y",
        "Y": "Y",
        "QE": "Q",
        "ME": "M",
        "YE": "Y",
    }
    if head not in mapping:
        raise SpecificationError(
            _i18n.message("time.freq_invalid", freq=freqstr, allowed=list(FREQS)),
            stata_code=198,
            key="time.freq_invalid",
        )
    return mapping[head]  # type: ignore[return-value]


def _integers(
    arr: np.ndarray,
    name: str,
    units: np.ndarray | None,
    unit_name: str,
    positions: IntArray | None,
) -> IntArray:
    """Integer values of a numeric column, refusing decimals (Stata error 451)."""
    if arr.dtype.kind in "iu":
        return arr.astype(np.int64)
    if arr.dtype.kind == "b":
        return arr.astype(np.int64)
    if arr.dtype.kind != "f":
        raise TimeAxisError(
            _i18n.message("time.not_numeric", name=name, dtype=str(arr.dtype)),
            key="time.not_numeric",
        )
    bad = ~np.isfinite(arr) | (arr != np.floor(arr))
    if bad.any():
        row = int(np.flatnonzero(bad)[0])
        raise TimeAxisError(
            _i18n.message(
                "time.not_integer",
                name=name,
                example=arr[row],
                row=_row(row, positions),
                unit=_where(row, units, unit_name),
            ),
            stata_code=451,
            key="time.not_integer",
        )
    return arr.astype(np.int64)


def encode_time(
    values: pd.Series | np.ndarray | pd.Index,
    *,
    name: str = "time",
    delta: int | None = None,
    freq: str | None = None,
    format: str | None = None,  # noqa: A002 - Stata's option name
    units: np.ndarray | None = None,
    unit_name: str = "id",
    positions: IntArray | None = None,
) -> EncodedTime:
    """Encode a time column as ranks on a regular calendar.

    Parameters
    ----------
    values
        Non-missing time values: integers, ``datetime64``, ``pandas.Period``, or integers coded
        for a Stata ``%t`` format (then pass ``format``).
    name
        Name of the time variable, used in messages.
    delta
        Calendar units between two consecutive periods (default 1, as in Stata).
    freq
        Calendar unit of dates: ``"D"``, ``"W"``, ``"M"``, ``"Q"``, ``"H"`` (half-year) or
        ``"Y"``.  Inferred for ``datetime64`` when omitted (the coarsest frequency on which all
        dates are period starts, or all period ends); carried by ``pandas.Period`` values.
    format
        Stata format (``"%tm"``, ``"%tq"``...) that says how integer values are coded.
    units, unit_name
        Identifier of each value (and its name), used to name the unit in messages.
    positions
        Row of the user's frame that holds each value, when ``values`` is a subset of its rows
        (messages then give the row of the frame).

    Returns
    -------
    EncodedTime
        The axis and an ``int64`` rank per row, with rank 0 at the earliest period.

    Raises
    ------
    TimeAxisError
        If the values are not numeric, dates or periods, are decimals, or are not on the
        calendar grid for ``delta``.
    SpecificationError
        If ``delta``, ``freq`` or ``format`` is invalid.

    Warns
    -----
    PyxtpanelWarning
        If all values are a multiple of some step larger than ``delta``, or if no date falls on a
        weekend (business days).
    """
    step = _check_delta(delta)
    freq_arg = _check_freq(freq)
    arr = np.asarray(values)
    is_dates = isinstance(getattr(values, "dtype", None), pd.PeriodDtype) or arr.dtype.kind == "M"

    kind: Kind
    used_freq: Freq | None
    if is_dates:
        days, kind, carried = _dates_from(values)
        if np.isnat(days).any():
            row = int(np.flatnonzero(np.isnat(days))[0])
            raise TimeAxisError(
                _i18n.message(
                    "time.missing",
                    name=name,
                    n=int(np.isnat(days).sum()),
                    rows=[_row(row, positions)],
                ),
                stata_code=459,
                key="time.missing",
            )
        used_freq = freq_arg or carried or _infer_freq(days)
        if used_freq == "D" and arr.dtype.kind == "M":
            full = arr.astype("datetime64[ns]")
            if np.any(full != full.astype("datetime64[D]").astype("datetime64[ns]")):
                example = full[np.flatnonzero(full != full.astype("datetime64[D]"))[0]]
                raise TimeAxisError(
                    _i18n.message("time.intraday", example=str(example)), key="time.intraday"
                )
        ordinal = _ordinals_from_dates(days, used_freq)
        if used_freq == "D" and freq_arg is None and kind == "datetime":
            _warn_business_days(days)
    elif format is not None:
        used_freq = freq_from_format(format)
        kind = "stata"
        ordinal = _integers(arr, name, units, unit_name, positions) + STATA_OFFSET[used_freq]
    else:
        if freq_arg is not None:
            raise SpecificationError(
                _i18n.message("time.freq_without_dates", freq=freq_arg, name=name),
                stata_code=198,
                key="time.freq_without_dates",
            )
        used_freq, kind = None, "integer"
        ordinal = _integers(arr, name, units, unit_name, positions)

    ordinal_min = int(np.min(ordinal)) if ordinal.size else 0
    offset = ordinal - ordinal_min
    if step > 1:
        off_grid = offset % step != 0
        if off_grid.any():
            row = int(np.flatnonzero(off_grid)[0])
            raise TimeAxisError(
                _i18n.message(
                    "time.not_aligned",
                    example=_example(values, row),
                    row=_row(row, positions),
                    unit=_where(row, units, unit_name),
                    delta=step,
                    first=_example(values, int(np.argmin(ordinal))),
                ),
                stata_code=451,
                key="time.not_aligned",
            )
    elif np.unique(offset).size >= 3:  # two periods alone say nothing about the step
        common = reduce(gcd, np.unique(offset).tolist(), 0)
        if common > 1:
            warnings.warn(
                _i18n.message("time.suggest_delta", delta=common),
                PyxtpanelWarning,
                stacklevel=3,
            )
    rank = (offset // step).astype(np.int64)
    n_steps = int(np.max(rank)) + 1 if rank.size else 0
    axis = TimeAxis(
        kind=kind,
        freq=used_freq,
        delta=step,
        ordinal_min=ordinal_min,
        n_steps=n_steps,
        stata_format=format,
    )
    return EncodedTime(axis=axis, rank=rank)


def _example(values: pd.Series | np.ndarray | pd.Index | Sequence[Any], row: int) -> object:
    item = np.asarray(values)[row]
    return item


def _warn_business_days(days: np.ndarray) -> None:
    unique = np.unique(days)
    if unique.size < 30:
        return
    weekday = (unique.astype(np.int64) + 3) % 7  # 0 = Monday
    if not np.any(weekday >= 5):
        warnings.warn(_i18n.message("time.business_days"), PyxtpanelWarning, stacklevel=4)


def check_span(axis: TimeAxis, n_units: int) -> None:
    """Raise if the lookup keys ``unit * n_steps + rank`` would overflow an int64.

    Raises
    ------
    TimeAxisError
        If ``n_units * n_steps`` is too large.
    """
    if n_units * max(axis.n_steps, 1) >= MAX_KEY:
        raise TimeAxisError(
            _i18n.message("time.span_too_large", n_steps=axis.n_steps, n_units=n_units),
            key="time.span_too_large",
        )
