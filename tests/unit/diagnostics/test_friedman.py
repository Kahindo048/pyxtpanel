"""Friedman test of cross-sectional independence: Spearman oracle and common periods."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError
from pyxtpanel.diagnostics import (
    friedman,
)


def declare(wide: np.ndarray) -> xt.PanelData:
    """Declare a ``(periods, units)`` array (NaN = absent) as a panel."""
    periods, units = wide.shape
    frame = pd.DataFrame(
        {
            "id": np.tile(np.arange(1, units + 1), periods),
            "t": np.repeat(np.arange(1, periods + 1), units),
            "v": wide.ravel(),
        }
    )
    return xt.xtset(frame.dropna(subset=["v"]).reset_index(drop=True), "id", "t")


def factor_data(seed: int, periods: int = 12, units: int = 8, load: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.normal(size=(periods, units)) + load * rng.normal(size=(periods, 1))


def test_friedman_by_hand_when_all_units_rank_alike() -> None:
    # Identical ranks everywhere: every Spearman correlation is 1, r = 1, T = 5, N = 3:
    # (T - 1) ((N - 1) r + 1) = 4 * 3 = 12 on N - 1 = 2 degrees of freedom.
    column = np.array([[3.0], [1.0], [4.0], [1.5], [9.0]])
    wide = np.hstack([column, 2 * column + 1, column**3])
    result = friedman(declare(wide), "v")
    assert result.statistic == pytest.approx(12.0, rel=1e-13) and result.df == (2.0,)
    assert result.pvalue == pytest.approx(float(stats.chi2.sf(12.0, 2)), rel=1e-12)


def test_friedman_matches_scipy_spearman() -> None:
    wide = factor_data(2)
    rhos = [stats.spearmanr(wide[:, i], wide[:, j])[0] for i, j in combinations(range(8), 2)]
    expected = (12 - 1) * ((8 - 1) * float(np.mean(rhos)) + 1)
    assert friedman(declare(wide), "v").statistic == pytest.approx(expected, rel=1e-12)


def test_friedman_uses_the_common_periods_and_the_shortest_unit_like_xtcsd() -> None:
    wide = factor_data(3)
    wide[:2, 0] = np.nan  # unit 1 starts late
    wide[10:, 1] = np.nan  # unit 2 ends early
    common = wide[2:10]
    rhos = [stats.spearmanr(common[:, i], common[:, j])[0] for i, j in combinations(range(8), 2)]
    shortest = 10  # units 1 and 2 have 10 periods each
    expected = (shortest - 1) * ((8 - 1) * float(np.mean(rhos)) + 1)
    result = friedman(declare(wide), "v")
    assert result.statistic == pytest.approx(expected, rel=1e-12)
    assert "8 periods where every unit is observed" in result.summary("en")
    wide[:, 3] = np.nan
    wide[:6, 3] = factor_data(3)[:6, 3]
    wide[:6, 4] = np.nan  # no period left with every unit
    with pytest.raises(DataError, match="no period has all 8 units observed"):
        friedman(declare(wide), "v")
