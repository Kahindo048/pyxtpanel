"""Slope homogeneity test: a unit-by-unit definition as oracle, options and refusals."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import run_test, slope_homogeneity


def make_panel(
    heterogeneous: bool, *, n: int = 40, t: int = 20, seed: int = 0, drop: float = 0.0
) -> xt.PanelData:
    rng = np.random.default_rng(seed)
    unit = np.repeat(np.arange(1, n + 1), t)
    spread = 1.0 if heterogeneous else 0.0
    slope = (1.0 + 0.5 * spread * rng.normal(size=n))[unit - 1]
    x, z, w = rng.normal(size=(3, n * t))
    effect = rng.normal(size=n)[unit - 1]
    y = effect + slope * x + 0.5 * z + 0.4 * w + rng.normal(size=n * t)
    frame = pd.DataFrame({"id": unit, "t": np.tile(np.arange(1, t + 1), n), "y": y, "x": x})
    frame["z"], frame["w"] = z, w
    if drop:
        frame = frame[rng.random(len(frame)) > drop].reset_index(drop=True)
    return xt.xtset(frame, "id", "t")


def by_definition(
    frame: pd.DataFrame, regressors: list[str], partial: tuple[str, ...] = (), *, ar: bool = False
) -> tuple[float, float, float]:
    """Delta, adjusted Delta and S from the formulas of Pesaran and Yamagata, unit by unit."""
    columns = ["y", *regressors]
    cleaned = []
    for _, g in frame.groupby("id"):
        z = np.column_stack([np.ones(len(g)), *[g[v].to_numpy() for v in partial]])
        data = g[columns].to_numpy()
        cleaned.append(data - z @ np.linalg.lstsq(z, data, rcond=None)[0])
    kp = 0 if ar else 1 + len(partial)
    xs = np.vstack([c[:, 1:] for c in cleaned])
    ys = np.concatenate([c[:, 0] for c in cleaned])
    pooled = np.linalg.solve(xs.T @ xs, xs.T @ ys)
    k, n = len(regressors), len(cleaned)
    weight_sum, weight_xy, parts = np.zeros((k, k)), np.zeros(k), []
    for c in cleaned:
        xi, yi = c[:, 1:], c[:, 0]
        e = yi - xi @ pooled
        s2 = float(e @ e) / (len(yi) - kp)
        bi = np.linalg.solve(xi.T @ xi, xi.T @ yi)
        weight_sum += xi.T @ xi / s2
        weight_xy += xi.T @ xi @ bi / s2
        parts.append((xi.T @ xi / s2, bi))
    b_wfe = np.linalg.solve(weight_sum, weight_xy)
    s = sum((bi - b_wfe) @ m @ (bi - b_wfe) for m, bi in parts)
    t_bar = len(frame) / n
    delta = np.sqrt(n) * (s / n - k) / np.sqrt(2 * k)
    var = 2 * k * (t_bar - k - kp - 1) / (t_bar - kp + 1)
    return float(delta), float(np.sqrt(n) * (s / n - k) / np.sqrt(var)), float(s)


def test_the_statistics_are_the_unit_by_unit_formulas() -> None:
    panel = make_panel(True)
    result = slope_homogeneity(panel, "y", ["x", "z"])
    delta, adjusted, swamy = by_definition(panel.to_frame(), ["x", "z"])
    assert result.details is not None
    assert result.details.loc["delta", "statistic"] == pytest.approx(delta, rel=1e-10)
    assert result.details.loc["delta_adjusted", "statistic"] == pytest.approx(adjusted, rel=1e-10)
    assert result.details.loc["swamy", "statistic"] == pytest.approx(swamy, rel=1e-10)
    assert result.statistic == pytest.approx(adjusted, rel=1e-10)
    assert result.pvalue == pytest.approx(2 * float(stats.norm.sf(abs(adjusted))), rel=1e-9)
    assert result.details.loc["swamy", "pvalue"] == pytest.approx(
        float(stats.chi2.sf(swamy, 2 * 39)), rel=1e-9
    )


def test_the_standard_statistic_is_selected_with_adjusted_false() -> None:
    panel = make_panel(True)
    standard = slope_homogeneity(panel, "y", ["x", "z"], adjusted=False)
    both = slope_homogeneity(panel, "y", ["x", "z"])
    assert both.details is not None
    assert standard.statistic == both.details.loc["delta", "statistic"]
    assert standard.statistic_name == "Delta" and both.statistic_name == "Delta_adj"


def test_common_slopes_are_not_rejected_and_different_slopes_are() -> None:
    assert not slope_homogeneity(make_panel(False, seed=3), "y", ["x"]).rejects(0.01)
    assert slope_homogeneity(make_panel(True, seed=3), "y", ["x"]).rejects(0.001)


def test_partial_noconstant_ar_and_averages_change_the_statistic_as_documented() -> None:
    panel = make_panel(True)
    base = slope_homogeneity(panel, "y", ["x"])
    assert slope_homogeneity(panel, "y", ["x"], partial=["z", "w"]).statistic != base.statistic
    assert slope_homogeneity(panel, "y", ["x"], constant=False).statistic != base.statistic
    # ar changes the degrees-of-freedom corrections, so the standard statistic moves too.
    with_ar = slope_homogeneity(panel, "y", ["x"], ar=True)
    assert with_ar.details is not None and base.details is not None
    assert with_ar.details.loc["delta", "statistic"] != base.details.loc["delta", "statistic"]
    averaged = slope_homogeneity(panel, "y", ["x"], crosssectional=["x", "y"])
    assert averaged.options["crosssectional"] == ("x", "y")
    assert averaged.statistic != base.statistic


def test_partialled_variables_and_the_ar_correction_follow_the_definition() -> None:
    panel = make_panel(True)
    frame = panel.to_frame()
    partial = slope_homogeneity(panel, "y", ["x"], partial=["z", "w"])
    delta, adjusted, swamy = by_definition(frame, ["x"], ("z", "w"))
    assert partial.details is not None
    assert partial.details.loc["delta", "statistic"] == pytest.approx(delta, rel=1e-10)
    assert partial.statistic == pytest.approx(adjusted, rel=1e-10)
    assert partial.details.loc["swamy", "statistic"] == pytest.approx(swamy, rel=1e-10)
    with_ar = slope_homogeneity(panel, "y", ["x"], partial=["w"], ar=True)
    assert with_ar.statistic == pytest.approx(
        by_definition(frame, ["x"], ("w",), ar=True)[1], rel=1e-10
    )


def test_an_unbalanced_panel_is_accepted() -> None:
    result = slope_homogeneity(make_panel(True, drop=0.1), "y", ["x"])
    assert np.isfinite(result.statistic) and result.n_obs < 800


def test_refusals() -> None:
    panel = make_panel(True)
    with pytest.raises(SpecificationError, match="at least one regressor"):
        slope_homogeneity(panel, "y", [])
    with pytest.raises(SpecificationError, match="both a tested variable and a partialled one"):
        slope_homogeneity(panel, "y", ["x"], partial=["x"])
    with pytest.raises(DataError, match="more observations than parameters"):
        slope_homogeneity(make_panel(True, t=3), "y", ["x", "z"], partial=["w"])
    twin = panel.to_frame().assign(x2=lambda f: 2 * f["x"])
    with pytest.raises(NumericalError, match="collinear"):
        slope_homogeneity(xt.xtset(twin, "id", "t"), "y", ["x", "x2"])


def test_by_name_and_alias_through_the_three_routes_in_both_languages() -> None:
    panel = make_panel(True)
    direct = slope_homogeneity(panel, "y", ["x"])
    assert run_test("xthst", panel, "y", ["x"]).statistic == direct.statistic
    assert xt.slope_heterogeneity(panel, "y", ["x"], test="swamy").statistic == direct.statistic
    assert panel.slope_heterogeneity("y", ["x"], test="slope_homogeneity").pvalue == direct.pvalue
    assert "Pesaran-Yamagata test of slope homogeneity" in direct.summary("en")
    assert "Test d'homogénéité des pentes" in direct.summary("fr")
    with pytest.raises(SpecificationError, match="choose a test"):
        xt.slope_heterogeneity(panel, "y", ["x"])


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


def test_too_few_periods_for_the_adjusted_statistic_is_refused_not_divided_by_zero() -> None:
    # one slope and a constant: the variance of the adjusted statistic needs T > 3
    panel = xt.xtset(
        _constant_unit_frame(units=12, periods=3).assign(y=lambda f: f["x"] ** 2), "id", "t"
    )
    with pytest.raises(DataError, match="needs more periods") as refusal:
        slope_homogeneity(panel, "y", ["x"])
    assert refusal.value.key == "diag.slope_homogeneity.too_short"


def hac_zero_bandwidth_by_definition(frame: pd.DataFrame) -> float:
    """S_HAC for one regressor and no lag: the Blomquist-Westerlund formulas in plain loops."""
    cleaned = []
    for _, g in frame.groupby("id"):
        data = g[["y", "x"]].to_numpy()
        data = data - data.mean(axis=0)  # the constant partialled out
        cleaned.append(data)
    xs = np.concatenate([c[:, 1] for c in cleaned])
    ys = np.concatenate([c[:, 0] for c in cleaned])
    pooled = float(xs @ ys) / float(xs @ xs)
    low = up = 0.0
    parts = []
    for c in cleaned:
        xi, yi = c[:, 1], c[:, 0]
        t = len(xi)
        scores = xi * (yi - xi * pooled)
        v = float(scores @ scores) / t
        q = float(xi @ xi) / t
        qvq = q * q / v
        low += t * qvq
        up += q / v * float(xi @ yi)
        parts.append((t, qvq, float(xi @ yi) / float(xi @ xi)))
    beta = up / low
    return sum(t * qvq * (b - beta) ** 2 for t, qvq, b in parts)


def test_the_hac_statistic_without_lags_is_the_white_form_from_the_definition() -> None:
    panel = make_panel(True, seed=3)
    result = slope_homogeneity(panel, "y", ["x"], bandwidth=0)
    s = hac_zero_bandwidth_by_definition(panel.frame)
    n = 40
    assert result.details is not None
    assert result.details.loc["swamy", "statistic"] == pytest.approx(s, rel=1e-10)
    assert result.details.loc["delta", "statistic"] == pytest.approx(
        np.sqrt(n) * (s / n - 1) / np.sqrt(2), rel=1e-10
    )
    assert result.options["bandwidth"] == 0.0 and result.options["kernel"] == "bartlett"


def test_whitening_and_a_bandwidth_imply_the_hac_version() -> None:
    panel = make_panel(True, seed=4)
    explicit = slope_homogeneity(panel, "y", ["x", "z"], hac=True, whitening=True)
    implied = slope_homogeneity(panel, "y", ["x", "z"], whitening=True)
    assert implied.statistic == explicit.statistic and "kernel" in implied.options
    plain = slope_homogeneity(panel, "y", ["x", "z"])
    assert "kernel" not in plain.options and plain.statistic != implied.statistic
    fixed = slope_homogeneity(panel, "y", ["x", "z"], bandwidth=2)
    assert fixed.options["bandwidth"] == 2.0 and fixed.options["bandwidth_fixed"] is True


@pytest.mark.parametrize("kernel", ["bartlett", "qs", "truncated"])
def test_every_kernel_gives_a_finite_statistic_and_a_bandwidth(kernel: str) -> None:
    result = slope_homogeneity(make_panel(True), "y", ["x", "z"], hac=True, kernel=kernel)  # type: ignore[arg-type]
    assert np.isfinite(result.statistic) and result.options["bandwidth"] >= 0  # type: ignore[operator]


def test_the_hac_version_runs_on_an_unbalanced_panel_and_with_three_regressors() -> None:
    result = slope_homogeneity(make_panel(True, drop=0.1), "y", ["x", "z", "w"], hac=True)
    assert np.isfinite(result.statistic)
    whitened = slope_homogeneity(make_panel(True), "y", ["x", "z", "w"], whitening=True)
    assert np.isfinite(whitened.statistic)  # may differ from xthst's (documented)


def test_lags_of_the_cross_sectional_averages_drop_the_rows_that_lack_them() -> None:
    panel = make_panel(True)
    none = slope_homogeneity(panel, "y", ["x"], crosssectional=["x"])
    lagged = slope_homogeneity(panel, "y", ["x"], crosssectional=["x"], cr_lags=2)
    assert none.n_obs == 800 and lagged.n_obs == 800 - 2 * 40
    assert lagged.options["cr_lags"] == (2,)
    # one value per variable, the last one applies to the following ones (as xthst)
    two = slope_homogeneity(panel, "y", ["x"], crosssectional=["x", "z"], cr_lags=[1, 2])
    assert two.options["cr_lags"] == (1, 2) and two.n_obs == 800 - 2 * 40
    same = slope_homogeneity(panel, "y", ["x"], crosssectional=["x", "z"], cr_lags=[1])
    assert same.options["cr_lags"] == (1, 1) and same.n_obs == 800 - 40


def test_hac_option_refusals() -> None:
    panel = make_panel(True)
    with pytest.raises(SpecificationError, match="hac and ar cannot be combined") as refusal:
        slope_homogeneity(panel, "y", ["x"], hac=True, ar=True)
    assert refusal.value.key == "diag.slope_homogeneity.hac_ar"
    with pytest.raises(SpecificationError, match="kernel must be"):
        slope_homogeneity(panel, "y", ["x"], hac=True, kernel="parzen")  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="bandwidth must be"):
        slope_homogeneity(panel, "y", ["x"], bandwidth=-1)
    with pytest.raises(SpecificationError, match="cr_lags must be"):
        slope_homogeneity(panel, "y", ["x"], crosssectional=["x"], cr_lags=-1)
