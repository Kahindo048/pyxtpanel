"""Harris-Tzavalis test: hand-computed value, explicit regressions as oracle, refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.diagnostics import harris_tzavalis, run_test


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


def walks(seed: int, n: int = 30, t: int = 10, rho: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    e = rng.normal(size=(n, t))
    y = np.empty((n, t))
    y[:, 0] = e[:, 0]
    for k in range(1, t):
        y[:, k] = rho * y[:, k - 1] + e[:, k]
    return y + rng.normal(size=(n, 1))


def test_hand_computed_without_constant() -> None:
    # Units (1, 2, 4) and (2, 1, 1): sum y_{t-1} y_t = 10 + 3, sum y_{t-1}^2 = 5 + 5, rho = 1.3;
    # T = 3, N = 2: mean 1, variance 2 / (T (T - 1)) / N = 1 / 6.
    result = harris_tzavalis(
        declare(np.array([[1.0, 2.0, 4.0], [2.0, 1.0, 1.0]])), "y", deterministics="none"
    )
    assert result.details is not None
    assert result.details.loc["rho", "value"] == pytest.approx(1.3, rel=1e-14)
    assert result.statistic == pytest.approx(0.3 * math.sqrt(6.0), rel=1e-13)
    assert result.pvalue == pytest.approx(float(stats.norm.cdf(0.3 * math.sqrt(6.0))), rel=1e-12)
    assert result.df == () and result.n_obs == 6 and result.n_groups == 2


def test_constant_matches_the_within_regression() -> None:
    y = walks(0, rho=0.8)
    result = harris_tzavalis(declare(y), "y")
    dep, lag = y[:, 1:], y[:, :-1]
    dep = dep - dep.mean(axis=1, keepdims=True)
    lag = lag - lag.mean(axis=1, keepdims=True)
    rho = (dep * lag).sum() / (lag**2).sum()
    n, t = y.shape
    mean = 1 - 3 / (t + 1)
    var = 3 * (17 * t**2 - 20 * t + 17) / (5 * (t - 1) * (t + 1) ** 3) / n
    assert result.details is not None
    assert result.details.loc["rho", "value"] == pytest.approx(rho, rel=1e-13)
    assert result.statistic == pytest.approx((rho - mean) / math.sqrt(var), rel=1e-12)


def test_trend_matches_unit_by_unit_detrending() -> None:
    y = walks(1, rho=0.7) + 0.4 * np.arange(10)
    result = harris_tzavalis(declare(y), "y", deterministics="trend")
    time = np.arange(2.0, 11.0)
    dep, lag = [], []
    for row in y:
        for series, store in ((row[1:], dep), (row[:-1], lag)):
            slope, intercept = np.polyfit(time, series, 1)
            store.append(series - slope * time - intercept)
    d, ell = np.array(dep), np.array(lag)
    rho = (d * ell).sum() / (ell**2).sum()
    assert result.details is not None
    assert result.details.loc["rho", "value"] == pytest.approx(rho, rel=1e-10)
    t = 10
    assert result.details.loc["mean_rho", "value"] == pytest.approx(1 - 15 / (2 * (t + 2)))


def test_demean_and_the_small_sample_adjustment() -> None:
    y = walks(2)
    centred = harris_tzavalis(declare(y - y.mean(axis=0)), "y")
    assert harris_tzavalis(declare(y), "y", demean=True).statistic == pytest.approx(
        centred.statistic, rel=1e-12
    )
    plain, adjusted = (
        harris_tzavalis(declare(y), "y"),
        harris_tzavalis(declare(y), "y", small_sample_t=True),
    )
    assert plain.details is not None and adjusted.details is not None
    assert plain.details.loc["rho", "value"] == adjusted.details.loc["rho", "value"]
    assert adjusted.details.loc["mean_rho", "value"] == pytest.approx(1 - 3 / 10)  # T = 9 -> T + 1


def test_stationary_panels_reject_and_random_walks_do_not() -> None:
    assert harris_tzavalis(declare(walks(3, 80, 10, rho=0.4)), "y").rejects(0.001)
    assert not harris_tzavalis(declare(walks(3, 80, 10, rho=1.0)), "y").rejects(0.001)


def test_refusals() -> None:
    y = walks(4)
    frame = declare(y).frame
    with pytest.raises(SpecificationError, match="deterministics must be"):
        harris_tzavalis(declare(y), "y", deterministics="drift")  # type: ignore[arg-type]
    with pytest.raises(DataError, match="strongly balanced panel: unit 3 is observed 9 times"):
        harris_tzavalis(
            xt.xtset(
                frame.drop(index=frame.index[(frame["id"] == 3) & (frame["t"] == 4)]), "id", "t"
            ),
            "y",
        )
    with pytest.raises(DataError, match="at least 4 observations per unit"):
        harris_tzavalis(declare(y[:, :3]), "y", deterministics="trend")
    with pytest.raises(DataError, match="at least 2 units"):
        harris_tzavalis(declare(y[:1]), "y")


def test_by_name_in_both_languages() -> None:
    result = run_test("harris_tzavalis", declare(walks(5)), "y")
    assert "Harris-Tzavalis unit-root test" in result.summary("en")
    assert "racine unitaire de Harris-Tzavalis" in result.summary("fr")


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


def test_a_constant_series_is_refused_not_given_a_nan_or_a_numerical_warning() -> None:
    from pyxtpanel._errors import NumericalError

    panel = xt.xtset(_constant_unit_frame(every_unit=True), "id", "t")
    with pytest.raises(NumericalError, match="no variation") as refusal:
        harris_tzavalis(panel, "y")
    assert refusal.value.key == "diag.harris_tzavalis.no_variation"


def test_one_constant_unit_among_others_still_gives_the_pooled_statistic() -> None:
    result = harris_tzavalis(xt.xtset(_constant_unit_frame(), "id", "t"), "y")
    assert math.isfinite(result.statistic)
