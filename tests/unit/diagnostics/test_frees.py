"""Frees test of cross-sectional independence: table, normal approximation and refusals."""

from __future__ import annotations

import math
from itertools import combinations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError
from pyxtpanel.diagnostics import (
    frees,
)
from pyxtpanel.diagnostics._shared.tables import load_table


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


def test_frees_matches_the_formula_and_the_critical_values_stata_prints() -> None:
    wide = factor_data(4)
    rhos = [stats.spearmanr(wide[:, i], wide[:, j])[0] for i, j in combinations(range(8), 2)]
    expected = 8 * (float(np.mean(np.square(rhos))) - 1 / (12 - 1))
    result = frees(declare(wide), "v")
    assert result.statistic == pytest.approx(expected, rel=1e-12)
    assert math.isnan(result.pvalue) and result.df == ()
    assert result.details is not None
    d = result.details["value"]
    # printed by Stata 17 for T = 12: 0.2136, 0.2838, 0.4252
    assert (d["cv_10"], d["cv_5"], d["cv_1"]) == (0.213585942, 0.283763514, 0.425186549)
    assert "critical values of Frees' Q" in result.summary("en")


def test_frees_table_rows_and_the_normal_approximation_for_long_panels() -> None:
    table = load_table("frees_q.csv")
    assert sorted(table["t"].unique()) == list(range(4, 31)) and len(table) == 81
    wide = factor_data(5, periods=40, units=6)
    result = frees(declare(wide), "v")
    m = 40.0
    seq = math.sqrt(
        32 / 25 * (m + 2) ** 2 / ((m - 1) ** 3 * (m + 1) ** 2)
        + 4 / 25 * (5 * m + 6) ** 2 * (m - 3) / (m * (m - 1) ** 2 * (m + 1) ** 2)
    )
    expected = min(1.0, 2 * float(stats.norm.sf(result.statistic / seq)))
    assert result.pvalue == pytest.approx(expected, rel=1e-12)
    independent = frees(declare(np.random.default_rng(6).normal(size=(40, 6))), "v")
    assert independent.pvalue <= 1.0  # capped even when Q is negative


def test_frees_needs_four_common_periods() -> None:
    with pytest.raises(DataError, match="fewer than 4 do"):
        frees(declare(factor_data(7, periods=3, units=5)), "v")
