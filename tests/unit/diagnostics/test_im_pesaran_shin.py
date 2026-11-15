"""Im-Pesaran-Shin test: tables, hand-checked statistics, unbalanced panels and refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import im_pesaran_shin, run_test
from pyxtpanel.diagnostics.unitroot.im_pesaran_shin import _critical_values, _table1, _table3


def walks(seed: int, n: int = 40, t: int = 20, rho: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    e = rng.normal(size=(n, t))
    y = np.empty((n, t))
    y[:, 0] = e[:, 0]
    for k in range(1, t):
        y[:, k] = rho * y[:, k - 1] + e[:, k]
    return y + rng.normal(size=(n, 1))


def declare(y: np.ndarray) -> xt.PanelData:
    n, t = y.shape
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n + 1), t),
            "t": np.tile(np.arange(1, t + 1), n),
            "y": y.ravel(),
        }
    )
    return xt.xtset(frame, "id", "t")


def test_table_1_rows_and_interpolation() -> None:
    assert _table1(6, True) == pytest.approx(
        {"mu_t": -1.520, "var_t": 1.745, "mu_tilde": -1.125, "var_tilde": 0.497}
    )
    assert _table1(4, True)["mu_t"] == -1.520  # below the first row: the first row
    assert _table1(20, True)["var_tilde"] == pytest.approx(0.592)
    low, high, middle = _table1(10, True), _table1(15, True), _table1(12, True)
    for key in middle:  # T = 12 is 40% of the way from the T = 10 row to the T = 15 row
        assert middle[key] == pytest.approx(0.6 * low[key] + 0.4 * high[key], abs=1e-12)
    assert _table1(12, False) == high  # without interpolation: the row at or above
    assert _table1(5000, True)["mu_t"] == -1.533  # the limit row


def test_table_2_critical_values() -> None:
    assert _critical_values(40, 20.0, False) == (-1.82, -1.73, -1.69)
    assert _critical_values(5, 5.0, False) == (-3.79, -2.76, -2.38)
    assert _critical_values(500, 200.0, True) == (-2.36, -2.31, -2.28)  # last row and column
    assert _critical_values(40, 20.0, True)[1] < _critical_values(40, 20.0, False)[1]


def test_table_3_cells_gaps_and_the_odd_stata_cell() -> None:
    assert _table3(0, 10, False, True) == (-1.504, 1.069)
    assert _table3(1, 20, False, True) == (-1.516, 0.915)
    assert _table3(2, 40, True, True) == (-2.120, 0.702)
    assert math.isnan(_table3(5, 10, False, True)[0])  # not in the paper
    assert math.isnan(_table3(9, 30, False, True)[0]) and math.isnan(
        _table3(-1, 30, False, True)[0]
    )
    # Stata's cell (trend, 8 lags, T = 60) is -2.204, out of line with its neighbours
    assert _table3(8, 60, True, True)[0] == -2.204
    low, high = _table3(1, 50, False, True), _table3(1, 60, False, True)
    mid = _table3(1, 55, False, True)
    assert mid[0] == pytest.approx((low[0] + high[0]) / 2) and mid[1] == pytest.approx(
        (low[1] + high[1]) / 2
    )
    assert _table3(1, 55, False, False) == high  # without interpolation: the column at or above
    assert _table3(1, 500, False, True) == _table3(1, 100, False, True)
    assert _table3(1, 8, False, True) == _table3(1, 10, False, True)


def test_statistics_without_lags_follow_the_published_formulas() -> None:
    y = walks(0)
    result = im_pesaran_shin(declare(y), "y")
    assert result.details is not None
    d = result.details["value"]
    taus, tildes = [], []
    for row in y:
        dy, level = np.diff(row), row[:-1]
        x = np.column_stack([level, np.ones(level.size)])
        beta = np.linalg.lstsq(x, dy, rcond=None)[0]
        resid = dy - x @ beta
        rmse = math.sqrt(resid @ resid / (dy.size - 2))
        se = rmse * math.sqrt(np.linalg.inv(x.T @ x)[0, 0])
        taus.append(beta[0] / se)
        tildes.append(beta[0] / (se / rmse * np.std(dy, ddof=1)))
    assert d["t_bar"] == pytest.approx(np.mean(taus), rel=1e-10)
    assert d["t_tilde_bar"] == pytest.approx(np.mean(tildes), rel=1e-10)
    # T = 20 exactly: Table 1 gives mean -1.395 and variance 0.592 for t-tilde, and
    # mean -1.522 and variance 0.851 for t
    assert result.statistic == pytest.approx(
        math.sqrt(40) * (np.mean(tildes) + 1.395) / math.sqrt(0.592), rel=1e-10
    )
    assert d["z_t_bar"] == pytest.approx(math.sqrt(40) * (np.mean(taus) + 1.522) / math.sqrt(0.851))
    assert result.pvalue == pytest.approx(float(stats.norm.cdf(result.statistic)))
    assert (d["cv_1"], d["cv_5"], d["cv_10"]) == (-1.82, -1.73, -1.69)


def test_w_t_bar_uses_table_3_for_the_lags_chosen() -> None:
    y = walks(1)
    result = im_pesaran_shin(declare(y), "y", lags=1)
    assert result.details is not None
    d = result.details["value"]
    assert (d["mu_bar"], d["var_bar"]) == (-1.516, 0.915) and d["lags"] == 1.0
    assert result.statistic == pytest.approx(
        math.sqrt(40) * (d["t_bar"] - d["mu_bar"]) / math.sqrt(d["var_bar"]), rel=1e-12
    )
    assert result.statistic_name == "W-t-bar"
    selected = im_pesaran_shin(declare(y), "y", lags=3, lag_selection="bic")
    assert selected.details is not None and 0.0 <= selected.details["value"]["lags"] <= 3.0


def test_stationary_panels_reject_and_random_walks_do_not() -> None:
    assert im_pesaran_shin(declare(walks(2, 40, 30, rho=0.5)), "y").rejects(0.001)
    assert not im_pesaran_shin(declare(walks(2, 40, 30, rho=1.0)), "y").rejects(0.001)


def test_demean_equals_subtracting_the_period_means_by_hand() -> None:
    y = walks(3)
    assert im_pesaran_shin(declare(y), "y", demean=True).statistic == pytest.approx(
        im_pesaran_shin(declare(y - y.mean(axis=0)), "y").statistic, rel=1e-12
    )


def test_an_unbalanced_panel_has_no_critical_values_and_a_short_one_no_p_value() -> None:
    y = walks(4, 30, 12)
    frame = declare(y).frame
    cut = frame[~((frame["id"] <= 10) & (frame["t"] > 8))]  # ten units end early: T = 8
    result = im_pesaran_shin(xt.xtset(cut, "id", "t"), "y")
    assert result.details is not None and math.isnan(result.details["value"]["cv_1"])
    assert math.isnan(result.pvalue) and "at least 10 observations" in result.summary("en")
    longer = frame[~((frame["id"] <= 10) & (frame["t"] > 10))]
    assert math.isfinite(im_pesaran_shin(xt.xtset(longer, "id", "t"), "y").pvalue)


def test_refusals() -> None:
    panel = declare(walks(5))
    for name, value in (("lag_selection", "fpe"), ("deterministics", "none")):
        with pytest.raises(SpecificationError, match=f"{name} must be one of"):
            im_pesaran_shin(panel, "y", **{name: value})  # type: ignore[arg-type]
    for bad in (-1, 1.5, True):
        with pytest.raises(SpecificationError, match="lags must be None or"):
            im_pesaran_shin(panel, "y", lags=bad)  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="give the largest number in lags"):
        im_pesaran_shin(panel, "y", lag_selection="aic")
    frame = panel.frame
    holed = xt.xtset(
        frame.drop(index=frame.index[(frame["id"] == 3) & (frame["t"] == 7)]), "id", "t"
    )
    with pytest.raises(DataError, match="consecutive observations: unit 3 has a gap"):
        im_pesaran_shin(holed, "y")
    with pytest.raises(NumericalError, match="not tabulated for T = 14 with 5 lag"):
        im_pesaran_shin(declare(walks(5, 20, 14)), "y", lags=5)
    with pytest.raises(NumericalError, match="not tabulated for T = 30 with 9 lag"):
        im_pesaran_shin(declare(walks(5, 20, 30)), "y", lags=9)


def test_by_name_in_both_languages() -> None:
    result = run_test("im_pesaran_shin", declare(walks(6)), "y")
    assert "Im-Pesaran-Shin unit-root test" in result.summary("en")
    assert "racine unitaire d'Im-Pesaran-Shin" in result.summary("fr")
