"""Levin-Lin-Chu test: kernels, bandwidths, adjustment table, lag selection and refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.diagnostics import levin_lin_chu, run_test
from pyxtpanel.diagnostics._shared.adf import select_lags
from pyxtpanel.diagnostics.unitroot import levin_lin_chu as module
from pyxtpanel.diagnostics.unitroot.levin_lin_chu import (
    _adjustments,
    _long_run_variance,
    _newey_west_bandwidth,
)


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


def series(seed: int, n: int = 25, t: int = 30, rho: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    e = rng.normal(size=(n, t))
    y = np.empty((n, t))
    y[:, 0] = e[:, 0]
    for k in range(1, t):
        y[:, k] = rho * y[:, k - 1] + e[:, k]
    return y + rng.normal(size=(n, 1))


def quantities(y: np.ndarray, **kwargs: object) -> pd.Series:
    result = levin_lin_chu(declare(y), "y", **kwargs)  # type: ignore[arg-type]
    assert result.details is not None
    return result.details["value"]


def test_the_adjusted_statistic_is_rebuilt_from_its_parts() -> None:
    result = levin_lin_chu(declare(series(0)), "y")
    assert result.details is not None
    d = result.details["value"]
    equations = d["t_tilde"] * result.n_groups
    expected = (
        d["t_delta"] - equations * d["s_bar"] / d["var_e"] * d["se_delta"] * d["mu_adj"]
    ) / d["sigma_adj"]
    assert result.statistic == pytest.approx(expected, rel=1e-12)
    assert result.pvalue == pytest.approx(float(stats.norm.cdf(expected)), rel=1e-12)
    assert d["t_delta"] == pytest.approx(d["delta"] / d["se_delta"], rel=1e-12)


def test_stationary_panels_reject_and_random_walks_do_not() -> None:
    assert levin_lin_chu(declare(series(1, 40, 40, rho=0.5)), "y").rejects(0.001)
    assert not levin_lin_chu(declare(series(1, 40, 40, rho=1.0)), "y").rejects(0.001)


def test_location_and_scale_do_not_change_the_statistic_and_demean_matches_by_hand() -> None:
    y = series(2)
    base = levin_lin_chu(declare(y), "y").statistic
    assert levin_lin_chu(declare(7.0 * y + 3.0), "y").statistic == pytest.approx(base, rel=1e-9)
    centred = levin_lin_chu(declare(y - y.mean(axis=0)), "y").statistic
    assert levin_lin_chu(declare(y), "y", demean=True).statistic == pytest.approx(
        centred, rel=1e-12
    )


def test_the_three_models_use_their_own_adjustment_row() -> None:
    y = series(3, 20, 60)
    # 58 equations per unit: 80% of the way from the T = 50 row to the T = 60 row
    mu_none = quantities(y, deterministics="none")["mu_adj"]
    mu_constant = quantities(y, deterministics="constant")["mu_adj"]
    mu_trend = quantities(y, deterministics="trend")["mu_adj"]
    assert mu_none == pytest.approx(0.001)
    assert mu_constant == pytest.approx(0.2 * -0.531 + 0.8 * -0.527)
    assert mu_trend == pytest.approx(0.2 * -0.614 + 0.8 * -0.598)


def test_adjustment_table_rows_and_interpolation() -> None:
    assert _adjustments(20.0, 2, True) == (-0.554, 0.919)  # below the first row
    assert _adjustments(25.0, 2, True) == (-0.554, 0.919)
    assert _adjustments(30.0, 3, True) == (-0.674, 0.949)  # exactly on a row
    assert _adjustments(1000.0, 1, True) == (0.0, 1.0)  # past the last finite row
    mu, sigma = _adjustments(27.5, 2, True)
    assert mu == pytest.approx((-0.554 - 0.546) / 2)
    assert sigma == pytest.approx((0.919 + 0.889) / 2)
    assert _adjustments(27.5, 2, False) == (-0.546, 0.889)  # the row at or above


def test_kernel_weights_by_hand() -> None:
    w = np.array([1.0, 2.0, 3.0, 4.0])
    gamma1, gamma2 = 2 + 6 + 12, 3 + 8  # sums of w_t w_{t-1} and w_t w_{t-2}
    total = 30.0
    assert _long_run_variance(w, "bartlett", 1) == pytest.approx((total + 2 * 0.5 * gamma1) / 4)
    assert _long_run_variance(w, "bartlett", 0) == pytest.approx(total / 4)
    bartlett2 = total + 2 * (1 - 1 / 3) * gamma1 + 2 * (1 - 2 / 3) * gamma2
    assert _long_run_variance(w, "bartlett", 2) == pytest.approx(bartlett2 / 4)
    weight1 = 1 - 6 * (1 / 3) ** 2 + 6 * (1 / 3) ** 3
    weight2 = 2 * (1 - 2 / 3) ** 3
    parzen = total + 2 * weight1 * gamma1 + 2 * weight2 * gamma2
    assert _long_run_variance(w, "parzen", 2) == pytest.approx(parzen / 4)


def test_quadratic_spectral_estimate_is_finite() -> None:
    w = np.random.default_rng(0).normal(size=40)
    value = _long_run_variance(w, "quadraticspectral", 4.5)
    assert value > 0 and math.isfinite(value)


def test_newey_west_bandwidth_is_an_integer_for_bartlett_and_parzen() -> None:
    rng = np.random.default_rng(1)
    w = np.convolve(rng.normal(size=300), [1.0, 0.7, 0.3], mode="valid")
    bartlett, parzen = _newey_west_bandwidth(w, "bartlett"), _newey_west_bandwidth(w, "parzen")
    assert bartlett == int(bartlett) and parzen == int(parzen) and bartlett >= 0
    spectral = _newey_west_bandwidth(w, "quadraticspectral")
    assert spectral != int(spectral)


def test_lag_selection_picks_the_dynamics_of_the_differences() -> None:
    rng = np.random.default_rng(2)
    e = rng.normal(size=600)
    d = np.zeros(600)
    for t in range(2, 600):  # the first difference is an AR(2): two lags are needed
        d[t] = 0.5 * d[t - 1] - 0.3 * d[t - 2] + e[t]
    y = np.cumsum(d)
    for criterion in ("aic", "bic", "hqic"):
        assert select_lags(y, 5, criterion) in (1, 2, 3)  # type: ignore[arg-type]
    assert select_lags(rng.normal(size=300).cumsum(), 4, "bic") <= 1


def test_without_a_constant_the_lag_selection_still_includes_one_as_stata_does(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: list[bool] = []
    real = module.select_lags

    def spy(*args: object, **kwargs: object) -> int:
        seen.append(bool(kwargs["constant"]))
        return real(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(module, "select_lags", spy)
    levin_lin_chu(
        declare(series(4, 10, 20)), "y", deterministics="none", lags=2, lag_selection="aic"
    )
    assert seen and all(seen)


def test_refusals() -> None:
    panel = declare(series(5))
    for name, value in (
        ("deterministics", "drift"),
        ("lag_selection", "fpe"),
        ("kernel", "tukey"),
    ):
        with pytest.raises(SpecificationError, match=f"{name} must be one of"):
            levin_lin_chu(panel, "y", **{name: value})  # type: ignore[arg-type]
    for bad in (-1, 1.5, True):
        with pytest.raises(SpecificationError, match="lags must be a non-negative integer"):
            levin_lin_chu(panel, "y", lags=bad)  # type: ignore[arg-type]
    for bad in ("fixed", -2, 2.5):
        with pytest.raises(SpecificationError, match="bandwidth must be"):
            levin_lin_chu(panel, "y", bandwidth=bad)  # type: ignore[arg-type]
    levin_lin_chu(panel, "y", kernel="quadraticspectral", bandwidth=2.5)  # a real value is allowed
    frame = panel.frame
    cut = xt.xtset(frame.drop(index=frame.index[(frame["id"] == 2) & (frame["t"] == 5)]), "id", "t")
    with pytest.raises(DataError, match="strongly balanced panel: unit 2 is observed 29 times"):
        levin_lin_chu(cut, "y")
    with pytest.raises(DataError, match="at least 6 observations per unit"):
        levin_lin_chu(declare(series(5, 10, 5)), "y", lags=2)


def test_by_name_in_both_languages() -> None:
    result = run_test("levin_lin_chu", declare(series(6)), "y")
    assert "Levin-Lin-Chu unit-root test" in result.summary("en")
    assert "racine unitaire de Levin-Lin-Chu" in result.summary("fr")


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


def test_a_constant_unit_is_refused_and_named() -> None:
    from pyxtpanel._errors import NumericalError

    panel = xt.xtset(_constant_unit_frame(), "id", "t")
    with pytest.raises(NumericalError, match="unit 1 has no variation") as refusal:
        levin_lin_chu(panel, "y", lags=1)
    assert refusal.value.key == "diag.levin_lin_chu.constant_unit"
