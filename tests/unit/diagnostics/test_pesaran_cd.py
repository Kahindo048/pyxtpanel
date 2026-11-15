"""Pesaran CD test: hand-computed values, an independent pairwise oracle, refusals."""

from __future__ import annotations

import math
from itertools import combinations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import pesaran_cd, run_test


def declare(frame: pd.DataFrame) -> xt.PanelData:
    return xt.xtset(frame, "id", "t")


def long_frame(wide: np.ndarray) -> pd.DataFrame:
    """Turn a ``(periods, units)`` array (NaN = absent) into a long frame."""
    periods, units = wide.shape
    frame = pd.DataFrame(
        {
            "id": np.tile(np.arange(1, units + 1), periods),
            "t": np.repeat(np.arange(1, periods + 1), units),
            "v": wide.ravel(),
        }
    )
    return frame.dropna(subset=["v"]).reset_index(drop=True)


def oracle(wide: np.ndarray) -> float:
    """CD from the definition, one pair at a time over the periods both units are observed."""
    n = wide.shape[1]
    total = 0.0
    for i, j in combinations(range(n), 2):
        both = ~np.isnan(wide[:, i]) & ~np.isnan(wide[:, j])
        rho = np.corrcoef(wide[both, i], wide[both, j])[0, 1]
        total += math.sqrt(both.sum()) * rho
    return math.sqrt(2.0 / (n * (n - 1))) * total


def test_hand_computed_on_three_units() -> None:
    # Three units, four periods: v1 = (1, 2, 3, 4), v2 = (2, 4, 6, 8) (rho = 1), v3 = (4, 3, 2, 1)
    # (rho = -1 with both): sum of rho = 1 - 1 - 1 = -1, CD = sqrt(2 * 4 / 6) * (-1).
    wide = np.array([[1.0, 2.0, 4.0], [2.0, 4.0, 3.0], [3.0, 6.0, 2.0], [4.0, 8.0, 1.0]])
    result = pesaran_cd(declare(long_frame(wide)), "v")
    assert result.statistic == pytest.approx(-math.sqrt(2 * 4 / 6), rel=1e-13)
    assert result.pvalue == pytest.approx(2 * stats.norm.sf(math.sqrt(2 * 4 / 6)), rel=1e-12)
    assert result.n_obs == 12 and result.n_groups == 3 and result.df == ()
    assert result.details is not None
    assert result.details.loc["average_correlation", "value"] == pytest.approx(-1 / 3)
    assert result.details.loc["average_abs_correlation", "value"] == pytest.approx(1.0)


def test_a_balanced_panel_matches_the_oracle() -> None:
    wide = np.random.default_rng(0).normal(size=(9, 7))
    assert pesaran_cd(declare(long_frame(wide)), "v").statistic == pytest.approx(
        oracle(wide), rel=1e-12
    )


def test_an_unbalanced_panel_uses_the_common_periods_of_every_pair() -> None:
    rng = np.random.default_rng(1)
    wide = rng.normal(size=(12, 9))
    wide[:3, 0] = np.nan  # late entry
    wide[9:, 1] = np.nan  # early exit
    wide[5, 2] = wide[7, 4] = wide[6, 4] = np.nan  # holes
    result = pesaran_cd(declare(long_frame(wide)), "v")
    assert result.statistic == pytest.approx(oracle(wide), rel=1e-12)
    assert "common periods" in result.summary("en")


def test_the_statistic_does_not_change_with_the_location_and_scale_of_each_unit() -> None:
    rng = np.random.default_rng(2)
    wide = rng.normal(size=(10, 6))
    changed = wide * rng.uniform(0.5, 3, size=6) + rng.normal(scale=10, size=6)
    assert pesaran_cd(declare(long_frame(changed)), "v").statistic == pytest.approx(
        pesaran_cd(declare(long_frame(wide)), "v").statistic, rel=1e-10
    )


def test_a_common_factor_rejects_and_independent_series_do_not() -> None:
    rng = np.random.default_rng(3)
    independent = rng.normal(size=(30, 25))
    common = independent + 1.5 * rng.normal(size=(30, 1))
    assert not pesaran_cd(declare(long_frame(independent)), "v").rejects(0.001)
    assert pesaran_cd(declare(long_frame(common)), "v").rejects(0.001)


def test_a_missing_value_removes_the_row() -> None:
    wide = np.random.default_rng(4).normal(size=(8, 5))
    frame = long_frame(wide)
    holed = frame.copy()
    holed.loc[7, "v"] = np.nan
    assert pesaran_cd(declare(holed), "v").statistic == pytest.approx(
        pesaran_cd(declare(frame.drop(index=7)), "v").statistic, rel=1e-13
    )


def test_refusals() -> None:
    wide = np.random.default_rng(5).normal(size=(6, 4))
    panel = declare(long_frame(wide))
    with pytest.raises(SpecificationError, match="min_overlap must be at least 3"):
        pesaran_cd(panel, "v", min_overlap=2)
    with pytest.raises(DataError, match="at least 3 units"):
        pesaran_cd(declare(long_frame(wide[:, :2])), "v")
    sparse = wide.copy()
    sparse[:4, 0] = np.nan  # unit 1 is observed 2 periods, so it shares 2 with every unit
    with pytest.raises(DataError, match="units 1 and 2 have 2 periods in common"):
        pesaran_cd(declare(long_frame(sparse)), "v")


def test_a_constant_pair_is_refused() -> None:
    wide = np.random.default_rng(6).normal(size=(6, 4))
    wide[:, 1] = 5.0
    with pytest.raises(NumericalError, match="no variation over their common periods"):
        pesaran_cd(declare(long_frame(wide)), "v")


def test_by_name_in_both_languages() -> None:
    wide = np.random.default_rng(7).normal(size=(8, 5))
    result = run_test("pesaran_cd", declare(long_frame(wide)), "v")
    assert "Pesaran CD test of cross-sectional dependence" in result.summary("en")
    assert "Test CD de Pesaran de dépendance transversale" in result.summary("fr")
