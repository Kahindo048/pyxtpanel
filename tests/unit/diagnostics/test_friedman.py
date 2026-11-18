"""Friedman test of cross-sectional independence: Spearman oracle and common periods."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, SpecificationError
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


def test_the_degrees_of_freedom_of_xtcsd_and_of_the_distribution_of_the_statistic() -> None:
    # Same by-hand panel: N = 3 units, T = 5 periods, statistic 12.
    column = np.array([[3.0], [1.0], [4.0], [1.5], [9.0]])
    panel = declare(np.hstack([column, 2 * column + 1, column**3]))
    default = friedman(panel, "v")
    periods = friedman(panel, "v", dof="periods")
    assert default.options["dof"] == "stata" and default.df == (2.0,)
    assert periods.options["dof"] == "periods" and periods.df == (4.0,)
    assert periods.statistic == default.statistic
    assert periods.pvalue == pytest.approx(float(stats.chi2.sf(12.0, 4)), rel=1e-12)
    # N != T: the default warns that the p-value of xtcsd is not the one of the distribution
    assert "as xtcsd does (N - 1)" in default.summary("en")
    assert "comme xtcsd (N - 1)" in default.summary("fr")
    assert "as xtcsd does" not in periods.summary("en")


def test_no_remark_when_n_equals_t_and_the_two_rules_then_agree() -> None:
    wide = factor_data(3, periods=9, units=9)
    default, periods = friedman(declare(wide), "v"), friedman(declare(wide), "v", dof="periods")
    assert default.pvalue == periods.pvalue and "as xtcsd does" not in default.summary("en")


def test_under_independence_the_statistic_is_a_chi_square_with_t_minus_1_degrees_of_freedom() -> (
    None
):
    # 6 units, 20 periods: the mean is T - 1 = 19 (not N - 1 = 5), and the p-values of the right
    # distribution are uniform while those of xtcsd are almost all zero.
    rng = np.random.default_rng(0)
    values, p_periods, p_stata = [], [], []
    for _ in range(300):
        panel = declare(rng.normal(size=(20, 6)))
        values.append(friedman(panel, "v").statistic)
        p_periods.append(friedman(panel, "v", dof="periods").pvalue)
        p_stata.append(friedman(panel, "v").pvalue)
    assert np.mean(values) == pytest.approx(19.0, rel=0.05)
    assert 0.02 <= np.mean(np.array(p_periods) < 0.05) <= 0.09
    assert np.mean(np.array(p_stata) < 0.05) > 0.9


def test_a_bad_dof_is_refused() -> None:
    with pytest.raises(SpecificationError, match="dof must be 'stata' or 'periods'"):
        friedman(declare(factor_data(4)), "v", dof="n")  # type: ignore[arg-type]


def _constant_unit_frame(
    units: int = 8, periods: int = 15, *, every_unit: bool = False
) -> pd.DataFrame:
    rng = np.random.default_rng(3)
    y = rng.normal(size=(units, periods)).cumsum(axis=1)
    y[: units if every_unit else 1] = 5.0
    return pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, units + 1), periods),
            "t": np.tile(np.arange(1, periods + 1), units),
            "y": y.ravel(),
            "x": rng.normal(size=units * periods),
        }
    )


def test_a_constant_unit_is_refused_because_its_ranks_are_all_tied() -> None:
    from pyxtpanel._errors import NumericalError
    from pyxtpanel.diagnostics import run_test

    panel = xt.xtset(_constant_unit_frame(), "id", "t")
    with pytest.raises(NumericalError, match="no variation") as refusal:
        run_test("friedman", panel, "y")
    assert refusal.value.key == "diag.pairwise.constant"
