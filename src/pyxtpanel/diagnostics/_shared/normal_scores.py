"""Normal scores of a sample, as the Shapiro tests of Stata compute them.

``swilk`` and ``sfrancia`` rank the data (ties receive the average rank, ``genrank``) and turn
the ranks into the approximate expected normal order statistics
``invnorm((rank - 0.375) / (n + 0.25))``; the statistic is a squared correlation between the data
and those scores or between the data and coefficients derived from them.  Without ties (option
``noties`` of Stata) the ranks are simply ``1..n`` in sorted order.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import stats

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray


def sorted_scores(x: FloatArray, *, ties: bool = True) -> tuple[FloatArray, FloatArray]:
    """Return the sorted data and the normal score of each value.

    Parameters
    ----------
    x
        ``(n,)`` finite values.
    ties
        Average rank for tied values (``True``, the default of Stata); ``False`` ranks tied
        values ``1..n`` in sorted order.
    """
    ordered = np.sort(x)
    n = ordered.size
    ranks = stats.rankdata(ordered, method="average") if ties else np.arange(1.0, n + 1.0)
    return ordered, np.asarray(stats.norm.ppf((ranks - 0.375) / (n + 0.25)), dtype=np.float64)


def squared_correlation(a: FloatArray, b: FloatArray) -> float:
    """Return the squared Pearson correlation of two vectors with variation."""
    ca = a - a.mean()
    cb = b - b.mean()
    return float((ca @ cb) ** 2 / ((ca @ ca) * (cb @ cb)))
