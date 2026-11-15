"""Breusch-Pagan LM test of cross-sectional independence: oracles and refusals."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.diagnostics import (
    lm_cross_section,
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


def test_lm_by_hand_on_three_perfectly_correlated_units() -> None:
    # v1 = (1, 2, 3, 4), v2 = 2 v1 (rho 1), v3 = (4, 3, 2, 1) (rho -1 with both): squared
    # correlations are all 1, so LM = 4 * 3 = 12 on 3 degrees of freedom.
    wide = np.array([[1.0, 2.0, 4.0], [2.0, 4.0, 3.0], [3.0, 6.0, 2.0], [4.0, 8.0, 1.0]])
    result = lm_cross_section(declare(wide), "v")
    assert result.statistic == pytest.approx(12.0, rel=1e-13) and result.df == (3.0,)
    assert result.pvalue == pytest.approx(float(stats.chi2.sf(12.0, 3)), rel=1e-12)


def test_lm_matches_a_pairwise_loop_on_an_unbalanced_panel() -> None:
    wide = factor_data(0)
    wide[:3, 0] = np.nan
    wide[9:, 1] = np.nan
    wide[5, 2] = wide[7, 4] = np.nan
    total = 0.0
    for i, j in combinations(range(wide.shape[1]), 2):
        both = ~np.isnan(wide[:, i]) & ~np.isnan(wide[:, j])
        total += both.sum() * np.corrcoef(wide[both, i], wide[both, j])[0, 1] ** 2
    result = lm_cross_section(declare(wide), "v")
    assert result.statistic == pytest.approx(total, rel=1e-12)
    assert "common periods" in result.summary("en")


def test_lm_warns_when_there_are_more_units_than_periods_and_refuses_bad_options() -> None:
    wide = factor_data(1, periods=6, units=9)
    assert "over-rejects" in lm_cross_section(declare(wide), "v").summary("en")
    with pytest.raises(SpecificationError, match="min_overlap must be at least 3"):
        lm_cross_section(declare(wide), "v", min_overlap=2)
    with pytest.raises(DataError, match="at least 3 units"):
        lm_cross_section(declare(wide[:, :2]), "v")
