"""Wooldridge autocorrelation test: an independent grid oracle, calendar rules and refusals."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import run_test, wooldridge


def make_frame(
    seed: int, units: int = 40, periods: int = 8, rho: float = 0.0, holes: float = 0.0
) -> pd.DataFrame:
    """A panel ``y = effect + 0.8 x - 0.4 z + u`` with AR(1) errors; cells may be removed."""
    rng = np.random.default_rng(seed)
    u = np.empty((units, periods))
    u[:, 0] = rng.normal(size=units)
    for t in range(1, periods):
        u[:, t] = rho * u[:, t - 1] + np.sqrt(1 - rho**2) * rng.normal(size=units)
    x, z = rng.normal(size=(2, units, periods))
    effect = rng.normal(size=(units, 1))
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, units + 1), periods),
            "t": np.tile(np.arange(1, periods + 1), units),
            "y": (effect + 0.8 * x - 0.4 * z + u).ravel(),
            "x": x.ravel(),
            "z": z.ravel(),
        }
    )
    if holes:
        frame = frame[rng.random(len(frame)) > holes].reset_index(drop=True)
    return frame


def oracle(frame: pd.DataFrame, names: list[str]) -> tuple[float, int, float]:
    """F, df_r and the slope from the full unit-by-period grid (NaN where nothing is observed)."""
    cube = np.stack(
        [frame.pivot(index="t", columns="id", values=v).to_numpy() for v in names], axis=-1
    )
    cube[np.isnan(cube).any(axis=-1)] = np.nan
    diff = cube[1:] - cube[:-1]
    ok = ~np.isnan(diff).any(axis=-1)
    dy, dx = diff[..., 0][ok], diff[..., 1:][ok]
    coef = np.linalg.lstsq(dx, dy, rcond=None)[0]
    res = np.full(ok.shape, np.nan)
    res[ok] = dy - dx @ coef
    y1, x1 = res[1:], res[:-1]
    both = ~np.isnan(y1) & ~np.isnan(x1)
    y1, x1 = np.where(both, y1, 0.0), np.where(both, x1, 0.0)
    slope = (x1 * y1).sum() / (x1 * x1).sum()
    score = (x1 * (y1 - slope * x1)).sum(axis=0)
    groups = int(both.any(axis=0).sum())
    variance = groups / (groups - 1) * (score**2).sum() / (x1 * x1).sum() ** 2
    return (slope + 0.5) ** 2 / variance, groups - 1, float(slope)


def test_matches_an_independent_grid_computation_on_a_balanced_panel() -> None:
    frame = make_frame(0, rho=0.4)
    result = wooldridge(xt.xtset(frame, "id", "t"), "y", ["x", "z"])
    f, df_r, slope = oracle(frame, ["y", "x", "z"])
    assert result.statistic == pytest.approx(f, rel=1e-10)
    assert result.df == (1.0, float(df_r))
    assert result.pvalue == pytest.approx(float(stats.f.sf(f, 1, df_r)), rel=1e-9)
    assert result.details is not None
    assert result.details.loc["residual_autocorrelation", "value"] == pytest.approx(slope)


def test_holes_make_the_differences_on_both_sides_missing() -> None:
    frame = make_frame(1, units=50, periods=10, rho=0.3, holes=0.15)
    result = wooldridge(xt.xtset(frame, "id", "t"), "y", ["x"])
    f, df_r, _ = oracle(frame, ["y", "x"])
    assert result.statistic == pytest.approx(f, rel=1e-10) and result.df == (1.0, float(df_r))


def test_a_missing_regressor_removes_the_row_like_a_hole() -> None:
    frame = make_frame(2, rho=0.2)
    with_nan = frame.copy()
    with_nan.loc[[10, 11, 50], "z"] = np.nan
    dropped = frame.drop(index=[10, 11, 50]).reset_index(drop=True)
    a = wooldridge(xt.xtset(with_nan, "id", "t"), "y", ["x", "z"])
    b = wooldridge(xt.xtset(dropped, "id", "t"), "y", ["x", "z"])
    assert a.statistic == pytest.approx(b.statistic, rel=1e-13) and a.df == b.df


def test_serial_correlation_is_detected_and_independence_is_not_rejected() -> None:
    correlated = wooldridge(xt.xtset(make_frame(3, 100, 8, rho=0.6), "id", "t"), "y", ["x", "z"])
    assert correlated.rejects(0.001)
    independent = wooldridge(xt.xtset(make_frame(3, 100, 8, rho=0.0), "id", "t"), "y", ["x", "z"])
    assert not independent.rejects(0.001)
    assert independent.details is not None
    assert abs(independent.details.loc["residual_autocorrelation", "value"] + 0.5) < 0.2


def test_unit_effects_and_shifts_of_the_variables_are_removed_by_differencing() -> None:
    frame = make_frame(4, rho=0.3)
    shifted = frame.copy()
    shifted["y"] += 5 * shifted["id"] + 100
    shifted["x"] += 3 * shifted["id"]
    a = wooldridge(xt.xtset(frame, "id", "t"), "y", ["x", "z"])
    b = wooldridge(xt.xtset(shifted, "id", "t"), "y", ["x", "z"])
    assert b.statistic == pytest.approx(a.statistic, rel=1e-8)


def test_a_collinear_regressor_does_not_change_the_residuals() -> None:
    frame = make_frame(5, rho=0.3)
    frame["x2"] = 2.0 * frame["x"]
    panel = xt.xtset(frame, "id", "t")
    assert wooldridge(panel, "y", ["x", "x2"]).statistic == pytest.approx(
        wooldridge(panel, "y", ["x"]).statistic, rel=1e-9
    )


def test_units_without_two_consecutive_differences_are_left_out_with_a_note() -> None:
    frame = make_frame(6, units=10, periods=6)
    frame = frame[~((frame["id"] == 1) & (frame["t"] > 2))]  # unit 1: one difference only
    result = wooldridge(xt.xtset(frame, "id", "t"), "y", ["x"])
    assert result.n_groups == 9 and result.df == (1.0, 8.0)
    assert "1 unit(s) lack two consecutive first differences" in result.summary("en")


def test_refusals() -> None:
    panel = xt.xtset(make_frame(7), "id", "t")
    with pytest.raises(SpecificationError, match="at least one regressor"):
        wooldridge(panel, "y", [])
    short = xt.xtset(make_frame(7, units=10, periods=2), "id", "t")  # one difference per unit
    with pytest.raises(NumericalError, match="at least 2"):
        wooldridge(short, "y", ["x"])
    flat = make_frame(7)
    flat["y"] = 1.0
    flat["x"] = 1.0  # all differences are zero
    with pytest.raises(NumericalError, match="residuals are all zero"):
        wooldridge(xt.xtset(flat, "id", "t"), "y", ["x"])
    empty = make_frame(7)
    empty["y"] = np.nan
    with pytest.raises(DataError, match="no observation is left"):
        wooldridge(xt.xtset(empty, "id", "t"), "y", ["x"])


def test_by_name_in_both_languages() -> None:
    panel = xt.xtset(make_frame(8, rho=0.5), "id", "t")
    result = run_test("wooldridge", panel, "y", ["x"])
    assert "Wooldridge test for autocorrelation in panel data" in result.summary("en")
    assert "autocorrélation dans un panel" in result.summary("fr")


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
    f = _constant_unit_frame(units=10, periods=12)
    f["y"] = np.random.default_rng(5).normal(size=len(f)).cumsum()
    f["x"] = 2.0 * f["y"] + 1.0
    with pytest.raises(NumericalError, match="all zero") as refusal:
        wooldridge(xt.xtset(f, "id", "t"), "y", ["x"])
    assert refusal.value.key == "diag.wooldridge.no_variation"
