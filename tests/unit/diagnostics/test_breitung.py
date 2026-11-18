"""Breitung test: invariances the transformation guarantees, a direct implementation, refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import breitung, run_test


def declare(y: np.ndarray) -> xt.PanelData:
    """Declare a ``(units, periods)`` array as a balanced panel."""
    n, t = y.shape
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n + 1), t),
            "t": np.tile(np.arange(1, t + 1), n),
            "y": y.ravel(),
        }
    )
    return xt.xtset(frame, "id", "t")


def series(seed: int, n: int = 25, t: int = 30, rho: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    e = rng.normal(size=(n, t))
    y = np.empty((n, t))
    y[:, 0] = e[:, 0]
    for k in range(1, t):
        y[:, k] = rho * y[:, k - 1] + e[:, k]
    return y + rng.normal(size=(n, 1))


def direct_constant(y: np.ndarray) -> float:
    """Breitung's statistic for a constant and no lag, written out unit by unit."""
    num = den = 0.0
    for row in y:
        dif = np.diff(row)
        sig2 = float(dif @ dif) / (row.size - 2)
        lag = row[:-1] - row[0]
        num += float(lag @ dif) / sig2
        den += float(lag @ lag) / sig2
    return num / math.sqrt(den)


def test_the_constant_statistic_follows_its_definition() -> None:
    y = series(0)
    result = breitung(declare(y), "y")
    assert result.statistic == pytest.approx(direct_constant(y), rel=1e-12)
    assert result.pvalue == pytest.approx(float(stats.norm.cdf(result.statistic)), rel=1e-12)
    assert result.n_obs == y.size and result.n_groups == 25 and result.df == ()


def test_adding_unit_intercepts_changes_nothing_with_a_constant_and_trends_with_a_trend() -> None:
    y = series(1)
    shift = np.random.default_rng(5).normal(scale=3.0, size=(y.shape[0], 1))
    slope = np.random.default_rng(6).normal(size=(y.shape[0], 1)) * np.arange(y.shape[1])
    for lags in (0, 2):
        base = breitung(declare(y), "y", lags=lags).statistic
        assert breitung(declare(y + shift), "y", lags=lags).statistic == pytest.approx(
            base, rel=1e-9
        )
        base_trend = breitung(declare(y), "y", deterministics="trend", lags=lags).statistic
        moved = breitung(declare(y + shift + slope), "y", deterministics="trend", lags=lags)
        assert moved.statistic == pytest.approx(base_trend, rel=1e-8)


def test_the_statistic_does_not_depend_on_the_unit_of_measurement() -> None:
    y = series(2)
    for kwargs in ({}, {"lags": 1}, {"deterministics": "trend"}, {"robust": True}):
        a = breitung(declare(y), "y", **kwargs).statistic  # type: ignore[arg-type]
        b = breitung(declare(1000.0 * y), "y", **kwargs).statistic  # type: ignore[arg-type]
        assert a == pytest.approx(b, rel=1e-9)


def test_demean_equals_subtracting_the_period_means_by_hand() -> None:
    y = series(3)
    assert breitung(declare(y), "y", demean=True, lags=1).statistic == pytest.approx(
        breitung(declare(y - y.mean(axis=0)), "y", lags=1).statistic, rel=1e-12
    )


def test_the_robust_statistic_is_reported_with_the_ordinary_one() -> None:
    y = series(4, n=10, t=40)
    result = breitung(declare(y), "y", robust=True)
    assert result.details is not None
    assert result.statistic == result.details.loc["lambda_robust", "value"]
    assert result.statistic_name == "lambda*"
    plain = breitung(declare(y), "y")
    assert plain.details is not None and "lambda_robust" not in plain.details.index
    assert result.details.loc["lambda", "value"] == pytest.approx(plain.statistic, rel=1e-12)


def test_it_rejects_stationary_panels_and_not_random_walks() -> None:
    assert breitung(declare(series(5, rho=0.5)), "y").rejects(0.001)
    assert not breitung(declare(series(5, rho=1.0)), "y").rejects(0.001)


def test_refusals() -> None:
    y = series(6)
    with pytest.raises(SpecificationError, match="deterministics must be"):
        breitung(declare(y), "y", deterministics="drift")  # type: ignore[arg-type]
    for bad in (-1, 1.5, True):
        with pytest.raises(SpecificationError, match="lags must be a non-negative integer"):
            breitung(declare(y), "y", lags=bad)  # type: ignore[arg-type]
    with pytest.raises(DataError, match="more periods, net of lags, than panels"):
        breitung(declare(series(7, n=30, t=20)), "y", robust=True)
    with pytest.raises(DataError):
        breitung(declare(y[:, :3]), "y", lags=2)
    with pytest.raises(NumericalError, match="no residual variation"):
        breitung(declare(np.ones((6, 10))), "y")


def test_by_name_alias_and_the_family_function_in_both_languages() -> None:
    panel = declare(series(8))
    direct = breitung(panel, "y", lags=1)
    assert run_test("breitung", panel, "y", lags=1).statistic == direct.statistic
    assert xt.unit_root(panel, "y", test="breitung", lags=1).statistic == direct.statistic
    assert "Breitung unit-root test" in direct.summary("en")
    assert "Test de racine unitaire de Breitung" in direct.summary("fr")
