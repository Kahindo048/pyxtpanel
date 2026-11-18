"""Hadri test: a direct implementation as oracle, invariances, the right tail, refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel import option_context
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import hadri, run_test


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


def series(seed: int, n: int = 30, t: int = 20, rho: float = 0.4) -> np.ndarray:
    rng = np.random.default_rng(seed)
    e = rng.normal(size=(n, t))
    y = np.empty((n, t))
    y[:, 0] = e[:, 0]
    for k in range(1, t):
        y[:, k] = rho * y[:, k - 1] + e[:, k]
    return y + rng.normal(size=(n, 1))


def direct(y: np.ndarray, *, trend: bool = False, robust: bool = False) -> tuple[float, float]:
    """The LM statistic and z by the definition, one unit at a time with polynomial fits."""
    n, t = y.shape
    degree = 1 if trend else 0
    divisor = t - 1 - degree
    lm_sum, ee, own = 0.0, 0.0, 0.0
    for row in y:
        e = row - np.polyval(np.polyfit(np.arange(1.0, t + 1), row, degree), np.arange(1.0, t + 1))
        s2 = float(np.sum(np.cumsum(e) ** 2))
        ee += float(e @ e)
        own += s2 / (float(e @ e) / divisor)
        lm_sum += s2
    lm = own / (n * t**2) if robust else lm_sum / (n * t**2) / (ee / (n * divisor))
    mean, var = (1 / 15, 11 / 6300) if trend else (1 / 6, 1 / 45)
    return lm, math.sqrt(n) * (lm - mean) / math.sqrt(var)


@pytest.mark.parametrize("trend", [False, True])
@pytest.mark.parametrize("robust", [False, True])
def test_the_statistic_follows_the_definition(trend: bool, robust: bool) -> None:
    y = series(0)
    result = hadri(declare(y), "y", deterministics="trend" if trend else "constant", robust=robust)
    lm, z = direct(y, trend=trend, robust=robust)
    assert result.details is not None
    assert result.details.loc["lm", "value"] == pytest.approx(lm, rel=1e-12)
    assert result.statistic == pytest.approx(z, rel=1e-11)
    assert result.pvalue == pytest.approx(float(stats.norm.sf(z)), rel=1e-11)
    assert result.n_obs == y.size and result.n_groups == 30


def test_with_a_kernel_the_variance_is_the_mean_long_run_variance() -> None:
    y = series(1)
    n, t = y.shape
    e = y - y.mean(axis=1, keepdims=True)
    lrv = []
    for row in e:
        value = float(row @ row) / t
        for lag in (1, 2):
            value += 2 * (1 - lag / 3) * float(row[lag:] @ row[:-lag]) / t
        lrv.append(value)
    lm = float(np.sum(np.cumsum(e, axis=1) ** 2)) / (n * t**2) / float(np.mean(lrv))
    result = hadri(declare(y), "y", kernel="bartlett", lags=2)
    assert result.details is not None and result.details.loc["lm", "value"] == pytest.approx(
        lm, rel=1e-12
    )


def test_the_default_kernel_length_is_one_lag() -> None:
    y = series(2)
    assert (
        hadri(declare(y), "y", kernel="parzen").statistic
        == hadri(declare(y), "y", kernel="parzen", lags=1).statistic
    )


def test_unit_intercepts_and_trends_and_the_scale_change_nothing() -> None:
    y = series(3)
    shift = np.random.default_rng(5).normal(scale=3.0, size=(y.shape[0], 1))
    slope = np.random.default_rng(6).normal(size=(y.shape[0], 1)) * np.arange(y.shape[1])
    assert hadri(declare(y + shift), "y").statistic == pytest.approx(
        hadri(declare(y), "y").statistic, rel=1e-9
    )
    assert hadri(
        declare(y + shift + slope), "y", deterministics="trend"
    ).statistic == pytest.approx(hadri(declare(y), "y", deterministics="trend").statistic, rel=1e-8)
    assert hadri(declare(100 * y), "y", robust=True).statistic == pytest.approx(
        hadri(declare(y), "y", robust=True).statistic, rel=1e-9
    )


def test_demean_equals_subtracting_the_period_means_by_hand() -> None:
    y = series(4)
    assert hadri(declare(y), "y", demean=True).statistic == pytest.approx(
        hadri(declare(y - y.mean(axis=0)), "y").statistic, rel=1e-12
    )


def test_the_null_is_stationarity_so_it_rejects_random_walks_in_the_right_tail() -> None:
    walks = hadri(declare(series(5, rho=1.0)), "y")
    stationary = hadri(declare(series(5, rho=0.0)), "y")
    assert walks.rejects(0.001) and walks.pvalue == pytest.approx(stats.norm.sf(walks.statistic))
    assert not stationary.rejects(0.001)


def test_refusals() -> None:
    y = series(6)
    with pytest.raises(SpecificationError, match="deterministics must be"):
        hadri(declare(y), "y", deterministics="none")  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="kernel must be"):
        hadri(declare(y), "y", kernel="tukey")  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="cannot be combined"):
        hadri(declare(y), "y", robust=True, kernel="bartlett")
    for bad in (1, 1.5, -1, True):
        kernel = "bartlett"
        if bad == 1:  # lags without a kernel
            kernel = None  # type: ignore[assignment]
        with pytest.raises(SpecificationError, match="lags must be"):
            hadri(declare(y), "y", kernel=kernel, lags=bad)  # type: ignore[arg-type]
    assert math.isfinite(hadri(declare(y), "y", kernel="quadraticspectral", lags=1.5).statistic)
    with pytest.raises(DataError):
        hadri(declare(y[:, :2]), "y", deterministics="trend")
    with pytest.raises(NumericalError, match="residual variance is zero"):
        hadri(declare(np.ones((6, 10))), "y")
    flat = y.copy()
    flat[0] = 4.0
    with pytest.raises(NumericalError, match="residual variance of unit 1 is zero"):
        hadri(declare(flat), "y", robust=True)


def test_by_name_alias_and_the_family_function_in_both_languages() -> None:
    panel = declare(series(7))
    direct_result = hadri(panel, "y")
    assert run_test("hadri", panel, "y").statistic == direct_result.statistic
    assert xt.unit_root(panel, "y", test="hadri").statistic == direct_result.statistic
    assert "Hadri LM stationarity test" in direct_result.summary("en")
    assert "stationnarité" in direct_result.summary("fr")


def test_a_flat_unit_is_named_in_the_robust_error_in_both_languages() -> None:
    y = (
        series(3, n=6, t=12)
        if "series" in globals()
        else np.random.default_rng(3).normal(size=(6, 12))
    )
    y[4] = 2.0  # the fifth unit, whose identifier is 5, is constant
    with pytest.raises(NumericalError, match=r"variance of unit 5 is zero") as english:
        hadri(declare(y), "y", robust=True)
    assert english.value.key == "diag.hadri.degenerate_unit"
    with (
        option_context(lang="fr"),
        pytest.raises(NumericalError, match=r"variance résiduelle de l'unité 5 est nulle"),
    ):
        hadri(declare(y), "y", robust=True)
