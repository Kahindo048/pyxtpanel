"""Panel Durbin-Watson and LBI: a loop over the units as oracle, gaps, delta and refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import panel_durbin_watson, run_test


def make_panel(
    rho: float = 0.0, *, n: int = 60, t: int = 10, drop: float = 0.0, step: int = 1, seed: int = 0
) -> tuple[xt.PanelData, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    errors = np.empty((n, t))
    errors[:, 0] = rng.normal(size=n)
    for k in range(1, t):
        errors[:, k] = rho * errors[:, k - 1] + np.sqrt(1 - rho**2) * rng.normal(size=n)
    unit = np.repeat(np.arange(1, n + 1), t)
    x = rng.normal(size=n * t)
    effect = rng.normal(size=n)[unit - 1]
    frame = pd.DataFrame(
        {
            "id": unit,
            "t": np.tile(np.arange(t) * step + 1, n),
            "x": x,
            "y": effect + x + errors.ravel(),
        }
    )
    if drop:
        keep = rng.random(len(frame)) > drop
        frame = frame[keep].reset_index(drop=True)
    return xt.xtset(frame, "id", "t", delta=step), frame


def by_loops(frame: pd.DataFrame, step: int = 1) -> tuple[float, float]:
    """d1 and LBI from the definitions, one unit at a time."""
    parts = []
    for _, g in frame.sort_values(["id", "t"]).groupby("id"):
        yc, xc = g["y"] - g["y"].mean(), g["x"] - g["x"].mean()
        parts.append((g["t"].to_numpy(), yc.to_numpy(), xc.to_numpy()))
    xs = np.concatenate([p[2] for p in parts])
    ys = np.concatenate([p[1] for p in parts])
    slope = float(xs @ ys) / float(xs @ xs)
    residuals = [(t, y - slope * x) for t, y, x in parts]
    total = sum(float(e @ e) for _, e in residuals)
    d1 = d2 = d3 = d4 = 0.0
    for t, e in residuals:
        for k in range(1, len(e)):
            previous = e[k - 1] if t[k] - t[k - 1] == step else 0.0
            d1 += (e[k] - previous) ** 2
        for k in range(len(e) - 1):
            if t[k + 1] - t[k] != step:
                d2 += e[k] ** 2
        d3 += e[0] ** 2
        d4 += e[-1] ** 2
    return d1 / total, (d1 + d2 + d3 + d4) / total


@pytest.mark.parametrize("drop", [0.0, 0.15])
@pytest.mark.parametrize("step", [1, 3])
def test_both_statistics_follow_the_definition_with_and_without_holes(
    drop: float, step: int
) -> None:
    panel, frame = make_panel(0.4, drop=drop, step=step, seed=2)
    result = panel_durbin_watson(panel, "y", ["x"])
    d1, lbi = by_loops(frame, step)
    assert result.details is not None
    assert result.statistic == pytest.approx(d1, rel=1e-10)
    assert result.details.loc["durbin_watson", "value"] == pytest.approx(d1, rel=1e-10)
    assert result.details.loc["baltagi_wu_lbi", "value"] == pytest.approx(lbi, rel=1e-10)


def test_independent_errors_give_about_two_and_positive_autocorrelation_gives_less() -> None:
    independent = panel_durbin_watson(make_panel(0.0, n=200)[0], "y", ["x"])
    assert independent.statistic == pytest.approx(2.0, abs=0.15)
    correlated = panel_durbin_watson(make_panel(0.7, n=200)[0], "y", ["x"])
    assert correlated.statistic < 1.0 and correlated.statistic < independent.statistic


def test_for_a_balanced_panel_the_lbi_adds_the_first_and_last_residuals() -> None:
    _, frame = make_panel(0.3, n=30, t=6, seed=5)
    d1, lbi = by_loops(frame)
    assert lbi > d1  # d3 + d4 > 0
    panel = xt.xtset(frame, "id", "t")
    assert panel_durbin_watson(panel, "y", ["x"]).details is not None


def test_there_is_no_p_value_and_the_remark_points_to_the_published_bounds() -> None:
    result = panel_durbin_watson(make_panel(0.5)[0], "y", ["x"])
    assert math.isnan(result.pvalue) and not result.rejects()
    assert "Bhargava, Franzini and Narendranathan (1982)" in result.summary("en")
    assert "bornes de Bhargava" in result.summary("fr")


def test_refusals() -> None:
    panel, frame = make_panel()
    with pytest.raises(SpecificationError, match="at least one regressor"):
        panel_durbin_watson(panel, "y", [])
    exact = xt.xtset(frame.assign(y=lambda f: 2 * f["x"]), "id", "t")
    with pytest.raises(NumericalError, match="all zero"):
        panel_durbin_watson(exact, "y", ["x"])


def test_by_name_alias_and_the_family_function_in_both_languages() -> None:
    panel, _ = make_panel(0.3)
    direct = panel_durbin_watson(panel, "y", ["x"])
    assert run_test("xtregar", panel, "y", ["x"]).statistic == direct.statistic
    assert xt.serial_correlation(panel, "y", ["x"], test="lbi").statistic == direct.statistic
    assert "Panel Durbin-Watson and Baltagi-Wu LBI statistics" in direct.summary("en")
    assert "Statistiques de Durbin-Watson de panel" in direct.summary("fr")


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
        panel_durbin_watson(xt.xtset(f, "id", "t"), "y", ["x"])
    assert refusal.value.key == "diag.panel_durbin_watson.zero_residuals"
