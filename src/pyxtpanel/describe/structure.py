"""``xtdescribe``: structure of the panel and participation patterns (E.1)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._results import Result
from pyxtpanel._types import FloatArray, Lang
from pyxtpanel.describe._spec import DescribeSpec, description
from pyxtpanel.describe.patterns import (
    check_options,
    participation_patterns,
    pattern_string,
)

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

#: Percentiles of ``T_i`` Stata returns in ``r()``.
PERCENTILES = (1, 5, 10, 25, 50, 75, 90, 95, 99)

_i18n.register(
    "describe",
    {
        "describe.xd.unique": {
            "en": "({id}*{time} uniquely identifies each observation)",
            "fr": "({id}*{time} identifie de façon unique chaque observation)",
        },
        "describe.xd.span": {"en": "{n} periods", "fr": "{n} périodes"},
        "describe.xd.distribution": {
            "en": "Distribution of T_i:   min      5%     25%       50%       75%     95%     max",
            "fr": "Distribution de T_i :  min      5%     25%       50%       75%     95%     max",
        },
        "describe.xd.header": {
            "en": "     Freq.  Percent    Cum. |  Pattern",
            "fr": "     Eff.   Pourc.    Cum. |  Motif",
        },
        "describe.xd.other": {"en": "(other patterns)", "fr": "(autres motifs)"},
        "describe.xd.each_column": {
            "en": "*Each column represents {dtv} periods.",
            "fr": "*Chaque colonne représente {dtv} périodes.",
        },
    },
)


def stata_quantile(values: FloatArray, percent: float) -> float:
    """Percentile as ``summarize, detail`` computes it.

    The average of two values when ``n·p/100`` is an integer.
    """
    return float(np.quantile(values, percent / 100.0, method="averaged_inverted_cdf"))


@dataclass(frozen=True, kw_only=True)
class XtDescribe(Result):
    """Result of :func:`xtdescribe`.

    Attributes
    ----------
    id_name, time_name
        Names of the panel variable and the time variable.
    n_units
        Number of units (``n``).
    unit_first, unit_second, unit_last
        Smallest, second smallest and largest identifier (as text).
    n_periods
        Distinct periods observed (``T``).
    span
        ``(tmax − tmin) / delta + 1`` periods.
    t_first, t_second, t_last
        First, second and last period observed (labels).
    delta_text
        Step of the calendar in words (``1 unit``, ``1 quarter``).
    ti
        ``summarize T_i, detail``: ``N``, ``mean``, ``sd``, ``min``, ``p1`` … ``p99``, ``max``.
    patterns
        Participation patterns shown: ``freq``, ``percent``, ``cum_percent``, ``pattern``.
    n_patterns
        Distinct patterns (shown or not).
    other_freq
        Units whose pattern is not shown.
    footer
        One symbol per column: ``X`` where some pattern is observed, ``.`` otherwise.
    periods_per_column
        Periods represented by one column of a pattern.
    """

    id_name: str
    time_name: str
    n_units: int
    unit_first: str
    unit_second: str
    unit_last: str
    n_periods: int
    span: int
    t_first: str
    t_second: str
    t_last: str
    delta_text: str
    ti: pd.Series
    patterns: pd.DataFrame
    n_patterns: int
    other_freq: int
    footer: str
    periods_per_column: int

    def to_frame(self) -> pd.DataFrame:
        """Return the patterns shown as a DataFrame."""
        return self.patterns

    def summary(self, lang: Lang | None = None) -> str:
        """Return the report in the layout of Stata's ``xtdescribe``."""

        def text(key: str, **fields: object) -> str:
            return _i18n.message(key, lang, **fields)

        t = self.ti
        first = f"{self.id_name:>8}:  {self.unit_first}, {self.unit_second}, ..., {self.unit_last}"
        second = f"{self.time_name:>8}:  {self.t_first}, {self.t_second}, ..., {self.t_last}"
        lines = [
            f"{first:<61}n ={self.n_units:>11}",
            f"{second:<61}T ={self.n_periods:>11}",
            f"{'':10}Delta({self.time_name}) = {self.delta_text}",
            f"{'':10}Span({self.time_name})  = {text('describe.xd.span', n=self.span)}",
            f"{'':10}{text('describe.xd.unique', id=self.id_name, time=self.time_name)}",
            "",
            text("describe.xd.distribution"),
            " " * 20 + f"{t['min']:6.0f}  {t['p5']:6.0f}  {t['p25']:6.0f}    {t['p50']:6.0f}    "
            f"{t['p75']:6.0f}  {t['p95']:6.0f}  {t['max']:6.0f}",
        ]
        if not len(self.patterns) and not self.other_freq:
            return "\n".join(lines)
        star = "*" if self.periods_per_column != 1 else ""
        width = max(len(self.footer), len(f"Pattern{star}"))
        lines += ["", text("describe.xd.header").replace("Pattern", f"Pattern{star}")]
        lines.append(" " + "-" * 27 + "+" + "-" * (width + 2))
        for row in self.patterns.itertuples():
            lines.append(
                f"{row.freq:>9d}{row.percent:>10.2f}{row.cum_percent:>8.2f} |  {row.pattern}"
            )
        if self.other_freq:
            share = 100 * self.other_freq / self.n_units
            lines.append(
                f"{self.other_freq:>9d}{share:>10.2f}{100:>8.2f} | {text('describe.xd.other')}"
            )
        lines.append(" " + "-" * 27 + "+" + "-" * (width + 2))
        lines.append(f"{self.n_units:>9d}{100:>10.2f}{'':8} |  {self.footer}")
        if star:
            lines += [
                "",
                " " + "-" * 27 + "-" + "-" * (width + 2),
                " " + text("describe.xd.each_column", dtv=self.periods_per_column),
            ]
        return "\n".join(lines)


@description(
    "xtdescribe",
    DescribeSpec(
        stata="xtdescribe",
        summary="Panel structure: units, periods, distribution of T_i, participation patterns.",
    ),
)
def xtdescribe(
    panel: PanelBase,
    *,
    variables: Sequence[str] | None = None,
    patterns: int = 9,
    width: int = 100,
) -> XtDescribe:
    """Describe the structure of a panel as Stata's ``xtdescribe`` does.

    Parameters
    ----------
    panel
        The declared panel.
    variables
        Count only the observations where all these variables are observed (the participation
        of the variables, not of the rows).
    patterns
        Number of patterns shown (default 9); ``0`` shows none; the rest is grouped.
    width
        Maximum columns of a pattern; beyond it a column represents several periods.

    Returns
    -------
    XtDescribe
        Counts, the distribution of ``T_i`` and the patterns.

    Raises
    ------
    DataError
        If fewer than two distinct periods are observed (code 459).
    SpecificationError
        If ``patterns`` or ``width`` is invalid.
    """
    check_options(patterns, width)
    sub = panel if variables is None else panel.select(complete=list(variables))
    index = sub.index
    periods = np.unique(index.rank)
    labels = sub.axis.labels(periods[[0, min(1, periods.size - 1), -1]])
    units = index.unit_labels
    counts = index.counts.astype(np.float64)
    ti = {"N": float(counts.size), "mean": float(counts.mean())}
    ti["sd"] = float(counts.std(ddof=1)) if counts.size > 1 else float("nan")
    ti["min"] = float(np.min(counts))
    for p in PERCENTILES:
        ti[f"p{p}"] = stata_quantile(counts, p)
    ti["max"] = float(np.max(counts))
    frame = pd.DataFrame(columns=["freq", "percent", "cum_percent", "pattern"])
    n_total, other, footer, dtv = 0, 0, "", 1
    if patterns != 0:
        found = participation_patterns(index, width)
        n_total = len(found.freq)
        shown = n_total if patterns >= 1000 else min(patterns, n_total)
        cum = np.cumsum(found.freq)
        frame = pd.DataFrame(
            {
                "freq": found.freq[:shown],
                "percent": 100.0 * found.freq[:shown] / index.n_units,
                "cum_percent": 100.0 * cum[:shown] / index.n_units,
                "pattern": [pattern_string(r) for r in found.counts[:shown]],
            }
        )
        other = int(cum[-1] - (cum[shown - 1] if shown else 0))
        footer = "".join("X" if on else "." for on in found.footer)
        dtv = found.periods_per_column
    else:
        participation_patterns(index, width)  # still refuses a single period (code 459)
    return XtDescribe(
        id_name=sub.id_name,
        time_name=sub.time_name,
        n_units=index.n_units,
        unit_first=str(units[0]),
        unit_second=str(units[min(1, units.size - 1)]),
        unit_last=str(units[-1]),
        n_periods=int(periods.size),
        span=int(periods[-1] - periods[0]) + 1,
        t_first=labels[0],
        t_second=labels[1],
        t_last=labels[2],
        delta_text=sub.axis.delta_text,
        ti=pd.Series(ti),
        patterns=frame,
        n_patterns=n_total,
        other_freq=other,
        footer=footer,
        periods_per_column=dtv,
    )
