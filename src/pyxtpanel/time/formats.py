"""Stata ``%t`` formats and the printable label of a calendar period.

Stata stores a period as an integer counted from 1960 (days since 1jan1960, weeks since 1960w1,
months since 1960m1, quarters since 1960q1, half-years since 1960h1; for ``%ty`` the integer is the
year itself).  The calendar ordinals of :mod:`pyxtpanel.time.axis` count from 1970, so a Stata
integer is moved by a fixed offset.

Stata's weeks (``%tw``) are not calendar weeks: a year has 52 of them, week ``w`` starts on 1
January plus ``7 * (w - 1)`` days, and the 52nd is 8 or 9 days long.  They are the frequency
``"SW"``; ``"W"`` is the calendar week of seven days that starts on a Monday.
"""

from __future__ import annotations

import re
from typing import Literal

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError

Freq = Literal["D", "W", "SW", "M", "Q", "H", "Y"]

#: Offset to add to a Stata-coded integer to obtain the calendar ordinal (counted from 1970).
STATA_OFFSET: dict[str, int] = {
    "D": -3653,  # 1960-01-01 is 3653 days before 1970-01-01
    "SW": -520,  # 10 years * 52 weeks
    "M": -120,  # 10 years * 12 months
    "Q": -40,
    "H": -20,
    "Y": -1970,  # %ty holds the year itself
}

_FORMAT = re.compile(r"^%-?t([dwmqhy])(?:\S*)$", re.IGNORECASE)
_LETTER_TO_FREQ: dict[str, Freq] = {"d": "D", "w": "SW", "m": "M", "q": "Q", "h": "H", "y": "Y"}

_MONTHS = ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")


def freq_from_format(fmt: str) -> Freq:
    """Return the frequency of a Stata ``%t`` format.

    Raises
    ------
    SpecificationError
        If the format is not one of ``%td %tw %tm %tq %th %ty`` (``%tc`` and ``%tC``, which count
        milliseconds, are not supported).
    """
    match = _FORMAT.match(fmt.strip())
    letter = match.group(1).lower() if match else ""
    if letter not in _LETTER_TO_FREQ:
        raise SpecificationError(
            _i18n.message("time.format_unsupported", fmt=fmt),
            stata_code=198,
            key="time.format_unsupported",
        )
    return _LETTER_TO_FREQ[letter]


def label_period(ordinal: int, freq: Freq) -> str:
    """Label of a calendar period as Stata prints it: ``2020q1``, ``2020m3``, ``01jan2020``."""
    if freq == "Y":
        return str(1970 + ordinal)
    if freq in ("M", "Q", "H"):
        divisor = {"M": 12, "Q": 4, "H": 2}[freq]
        year, index = 1970 + ordinal // divisor, ordinal % divisor + 1
        return f"{year}{freq.lower()}{index}"
    if freq == "SW":
        return f"{1970 + ordinal // 52}w{ordinal % 52 + 1}"
    if freq == "W":
        day = np.datetime64(ordinal * 7 - 3, "D")
        return _day_label(day)
    return _day_label(np.datetime64(ordinal, "D"))


def _day_label(day: np.datetime64) -> str:
    text = str(day)  # YYYY-MM-DD
    year, month, number = text.split("-")
    return f"{number}{_MONTHS[int(month) - 1]}{year}"
