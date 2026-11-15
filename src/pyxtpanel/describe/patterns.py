"""Participation patterns of ``xtdescribe``: which units are observed in which periods.

Each column of a pattern is a period (or a group of ``dtv`` periods when there are more periods
than ``width``); its value is the number of observations of the unit in the column: ``.`` for 0,
the digit for 1 to 9, ``X`` for 10 or more.  Units with the same pattern are grouped.  Rules and
sort order are those of Stata's ``xtdescribe.ado`` (ANALYSE_STATA.md, section 2).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel._types import IntArray
from pyxtpanel.panel.index import PanelIndex

_i18n.register(
    "describe",
    {
        "describe.few_periods": {
            "en": "the time variable must have multiple distinct nonmissing values.",
            "fr": "la variable de temps doit avoir plusieurs valeurs distinctes non manquantes.",
        },
        "describe.bad_width": {
            "en": "width must be an integer greater than 2, got {width}.",
            "fr": "width doit être un entier supérieur à 2, reçu {width}.",
        },
        "describe.bad_patterns": {
            "en": "patterns must be a non-negative integer, got {patterns}.",
            "fr": "patterns doit être un entier positif ou nul, reçu {patterns}.",
        },
    },
)


@dataclass(frozen=True)
class Patterns:
    """Participation patterns of a panel.

    Attributes
    ----------
    counts
        ``(n_patterns, n_columns)`` observations per column, most frequent pattern first, ties
        broken by the counts in increasing order (Stata's order).
    freq
        Units with each pattern.
    footer
        Per column, whether any pattern is observed there.
    periods_per_column
        ``dtv``: periods represented by a column (1 unless the panel is wider than ``width``).
    """

    counts: IntArray
    freq: IntArray
    footer: np.ndarray
    periods_per_column: int


def symbol(count: int) -> str:
    """``.`` for 0, the digit for 1 to 9, ``X`` for 10 or more."""
    return "." if count == 0 else str(count) if count <= 9 else "X"


def pattern_string(row: np.ndarray) -> str:
    """Text of one pattern (one symbol per column)."""
    return "".join(symbol(int(c)) for c in row)


def check_options(patterns: int, width: int) -> None:
    """Validate the ``patterns`` and ``width`` options like Stata does."""
    if patterns < 0:
        raise SpecificationError(
            _i18n.message("describe.bad_patterns", patterns=patterns), key="describe.bad_patterns"
        )
    if patterns != 0 and width <= 2:
        raise SpecificationError(
            _i18n.message("describe.bad_width", width=width), key="describe.bad_width"
        )


def periods_per_column(distinct: IntArray, width: int) -> int:
    """``dtv`` of ``xtdescribe``: smallest observed step, widened to fit ``width`` columns."""
    # STATA-QUIRK: xtdescribe sizes a column from the smallest observed step, not from delta
    smallest = int(np.min(np.diff(distinct)))
    span = int(distinct[-1] - distinct[0])
    return span // (width - 1) if span // smallest + 1 > width else smallest


def participation_patterns(index: PanelIndex, width: int = 100) -> Patterns:
    """Compute the patterns of the sorted index.

    Raises
    ------
    DataError
        If fewer than two distinct periods are observed (Stata code 459).
    """
    distinct = np.unique(index.rank)
    if distinct.size < 2:
        raise DataError(
            _i18n.message("describe.few_periods"), stata_code=459, key="describe.few_periods"
        )
    dtv = periods_per_column(distinct, width)
    relative = index.rank - distinct[0]
    last = int(distinct[-1] - distinct[0])
    n_columns = int(np.ceil((last + 0.5) / dtv))
    column = np.ceil((relative + 0.5) / dtv).astype(np.int64) - 1
    table = np.bincount(
        index.unit * n_columns + column, minlength=index.n_units * n_columns
    ).reshape(index.n_units, n_columns)
    rows, freq = np.unique(table, axis=0, return_counts=True)
    # most frequent first; equal frequencies by the counts per column, increasing
    keys = [rows[:, j] for j in range(n_columns - 1, -1, -1)] + [-freq]
    order = np.lexsort(keys)
    rows, freq = rows[order], freq[order]
    return Patterns(
        counts=rows.astype(np.int64),
        freq=freq.astype(np.int64),
        footer=(rows > 0).any(axis=0),
        periods_per_column=dtv,
    )
