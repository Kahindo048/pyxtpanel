"""Kao tests: identities the statistics must satisfy, a direct implementation, refusals."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import kao, run_test
from pyxtpanel.diagnostics._shared.long_run import (
    kernel_weight,
    long_run_covariance,
    newey_west_lags,
)

DATA = Path(__file__).resolve().parents[2] / "golden" / "data"


def simulate(seed: int, *, n: int = 20, t: int = 40, cointegrated: bool = True) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    shocks = rng.normal(size=(n, t, 2))
    shocks[:, :, 1] += 0.5 * shocks[:, :, 0]
    x = np.cumsum(shocks[:, :, 1], axis=1)
    u = np.cumsum(shocks[:, :, 0], axis=1) if not cointegrated else shocks[:, :, 0] * 0.8
    y = rng.normal(scale=2.0, size=(n, 1)) + x + u
    return pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n + 1), t),
            "t": np.tile(np.arange(1, t + 1), n),
            "y": y.ravel(),
            "x": x.ravel(),
        }
    )


def declare(frame: pd.DataFrame) -> xt.PanelData:
    return xt.xtset(frame, "id", "t")


def test_it_rejects_a_cointegrated_panel_and_not_a_spurious_one() -> None:
    assert kao(declare(simulate(0)), "y", ["x"]).rejects(0.001)
    assert not kao(declare(simulate(0, cointegrated=False)), "y", ["x"]).rejects(0.001)


def test_the_five_statistics_and_the_choice_of_the_reported_one() -> None:
    panel = declare(simulate(1))
    result = kao(panel, "y", ["x"])
    assert result.details is not None
    assert list(result.details.index) == [
        "modified_df",
        "df",
        "adf",
        "unadjusted_modified_df",
        "unadjusted_df",
    ]
    assert result.statistic_name == "adf" and result.n_groups == 20 and result.n_obs == 20 * 38
    assert result.statistic == result.details.loc["adf", "statistic"]
    for name in result.details.index:
        chosen = kao(panel, "y", ["x"], statistic=name)  # type: ignore[arg-type]
        assert chosen.statistic == result.details.loc[name, "statistic"]
        assert chosen.pvalue == pytest.approx(result.details.loc[name, "pvalue"])


def test_the_unadjusted_statistics_follow_their_formulas() -> None:
    frame = simulate(2)
    result = kao(declare(frame), "y", ["x"])
    assert result.details is not None
    n, periods = 20, 40
    y = frame["y"].to_numpy().reshape(n, periods)[:, 1:]
    x = frame["x"].to_numpy().reshape(n, periods)[:, 1:]
    yc, xc = y - y.mean(axis=1, keepdims=True), x - x.mean(axis=1, keepdims=True)
    slope = float((xc * yc).sum() / (xc * xc).sum())
    e = yc - slope * xc
    lag, cur = e[:, :-1], e[:, 1:]
    rho = float((lag * cur).sum() / (lag * lag).sum())
    t = periods - 2
    resid = cur - rho * lag
    t_rho = (rho - 1.0) * math.sqrt(float((lag * lag).sum()) / (float((resid**2).sum()) / (n * t)))
    expected_rho = (math.sqrt(n) * t * (rho - 1.0) + 3.0 * math.sqrt(n)) / math.sqrt(51.0 / 5.0)
    expected_t = math.sqrt(5.0 / 4.0) * t_rho + math.sqrt(15.0 * n / 8.0)
    assert result.details.loc["unadjusted_modified_df", "statistic"] == pytest.approx(
        expected_rho, rel=1e-12
    )
    assert result.details.loc["unadjusted_df", "statistic"] == pytest.approx(expected_t, rel=1e-12)


def test_without_kernel_lags_the_adjusted_statistics_are_the_unadjusted_ones() -> None:
    # Sigma = Omega, so the variance ratio is 1 and the corrections of the formulas vanish.
    result = kao(declare(simulate(3)), "y", ["x"], kernel_lags=0)
    assert result.details is not None
    table = result.details["statistic"]
    assert table["modified_df"] == pytest.approx(table["unadjusted_modified_df"], rel=1e-12)
    assert table["df"] == pytest.approx(table["unadjusted_df"], rel=1e-12)


def test_the_scale_does_not_change_the_statistics_and_demean_is_the_period_mean_removed() -> None:
    # The innovations come from a pooled AR(1) without constant, so unit means do matter (as in
    # Stata); a common change of units does not.
    frame = simulate(4)
    base = kao(declare(frame), "y", ["x"]).details
    scaled = frame.copy()
    scaled[["y", "x"]] *= 1000.0
    moved = kao(declare(scaled), "y", ["x"]).details
    assert base is not None and moved is not None
    assert moved["statistic"].to_numpy() == pytest.approx(base["statistic"].to_numpy(), rel=1e-7)
    by_hand = frame.copy()
    for column in ("y", "x"):
        by_hand[column] -= by_hand.groupby("t")[column].transform("mean")
    demeaned = kao(declare(frame), "y", ["x"], demean=True).details
    direct = kao(declare(by_hand), "y", ["x"]).details
    assert demeaned is not None and direct is not None
    assert demeaned["statistic"].to_numpy() == pytest.approx(
        direct["statistic"].to_numpy(), rel=1e-10
    )


def test_the_automatic_lags_are_those_stata_reports_on_the_frozen_panels() -> None:
    # Stata 17: xtcointtest kao y x z on sim_coint_yes / sim_coint_no reports Lags 2.57 and 2.20.
    for name, lags in (("coint_yes", 77 / 30), ("coint_no", 66 / 30)):
        frame = pd.read_csv(DATA / f"{name}.csv")
        result = kao(declare(frame), "y", ["x", "z"])
        assert result.options["kernel_lags"] == pytest.approx(lags, rel=1e-12)


def test_the_newey_west_rule_and_the_kernels() -> None:
    rng = np.random.default_rng(5)
    white = rng.normal(size=(200, 2))
    lags, bandwidth = newey_west_lags(white, "bartlett")
    assert lags <= int(4 * (200 / 100) ** (2 / 9)) and bandwidth >= 0
    assert kernel_weight("bartlett", 1, 4.0) == pytest.approx(0.75)
    assert kernel_weight("parzen", 1, 4.0) == pytest.approx(1 - 6 / 16 + 6 / 64)
    assert kernel_weight("parzen", 3, 4.0) == pytest.approx(2 * (1 / 4) ** 3)
    assert kernel_weight("quadraticspectral", 0, 4.0) == 1.0
    fixed = long_run_covariance(white[None, :, :], "bartlett", 2)
    assert fixed is not None and fixed.lags == 2 and fixed.bandwidth == 3.0
    centred = white - white.mean(axis=0)
    assert fixed.sigma == pytest.approx(centred.T @ centred / 200, rel=1e-12)
    zero = long_run_covariance(white[None, :, :], "bartlett", 0)
    assert zero is not None and zero.omega == pytest.approx(zero.sigma, rel=1e-14)


def test_the_selection_of_augmented_lags_picks_a_number_in_range() -> None:
    panel = declare(simulate(6))
    for criterion in ("aic", "bic", "hqic"):
        result = kao(panel, "y", ["x"], adf_lags=4, adf_lag_selection=criterion)  # type: ignore[arg-type]
        assert 1 <= int(result.options["adf_lags"]) <= 4  # type: ignore[call-overload]


def test_refusals() -> None:
    panel = declare(simulate(7))
    with pytest.raises(SpecificationError, match="at least one regressor"):
        kao(panel, "y", [])
    for name, value in (("statistic", "z"), ("kernel", "tukey"), ("adf_lag_selection", "ic")):
        with pytest.raises(SpecificationError, match=f"{name} must be one of"):
            kao(panel, "y", ["x"], **{name: value})  # type: ignore[arg-type]
    for name, value in (
        ("adf_lags", 0),
        ("kernel_lags", -1),
        ("kernel_lags", 1.5),
        ("adf_lags", True),
    ):
        with pytest.raises(SpecificationError, match=f"{name} must be a non-negative integer"):
            kao(panel, "y", ["x"], **{name: value})  # type: ignore[arg-type]
    frame = simulate(7)
    with pytest.raises(DataError):
        kao(declare(frame[frame["t"] <= 5]), "y", ["x"], adf_lags=3)
    frame["twin"] = 2.0 * frame["x"]
    with pytest.raises(NumericalError, match="no usable variation"):
        kao(declare(frame), "y", ["x", "twin"])


def test_by_name_the_family_function_and_the_method_in_both_languages() -> None:
    panel = declare(simulate(8))
    direct = kao(panel, "y", ["x"])
    assert run_test("kao", panel, "y", ["x"]).statistic == direct.statistic
    assert xt.cointegration(panel, "y", ["x"], test="kao").statistic == direct.statistic
    assert panel.cointegration("y", ["x"], test="kao").pvalue == direct.pvalue
    with pytest.raises(SpecificationError, match="choose a test"):
        xt.cointegration(panel, "y", ["x"])
    assert "Kao test for cointegration" in direct.summary("en")
    assert "Test de cointégration de Kao" in direct.summary("fr")


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


def test_a_regressor_that_fits_the_dependent_variable_perfectly_is_refused() -> None:
    f = _constant_unit_frame(units=10, periods=20)
    f["y"] = np.random.default_rng(5).normal(size=len(f)).cumsum()
    f["x"] = 2.0 * f["y"] + 1.0
    with pytest.raises(NumericalError, match="no usable variation") as refusal:
        kao(xt.xtset(f, "id", "t"), "y", ["x"], kernel_lags=2)
    assert refusal.value.key == "diag.kao.degenerate"
