"""CIPS test: R's plm values, a direct implementation, the tables, the options, refusals."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import cips, run_test

DATA = Path(__file__).resolve().parents[2] / "golden" / "data"

#: plm 2.6.7, R 4.6.1: cipstest(pdata$y, lags, type, truncated = FALSE) on the frozen balanced
#: panels: (panel, lags, deterministics) -> CIPS statistic (the mean of the unit t ratios).  plm
#: does not truncate by default either.
PLM = {
    ("cips_factor", 1, "constant"): -2.251894190806059,
    ("cips_factor", 2, "trend"): -2.731517699316683,
    ("ur_mixed", 1, "constant"): -1.752628084775661,
    ("ur_long", 3, "trend"): -2.78253615191444,
}


def declare(frame: pd.DataFrame) -> xt.PanelData:
    return xt.xtset(frame, "id", "t")


def simulate(
    seed: int, *, n: int = 20, t: int = 30, rho: float = 1.0, factor: float = 1.0
) -> pd.DataFrame:
    """Series with an innovation factor common to the units and a root ``rho``."""
    rng = np.random.default_rng(seed)
    common = rng.normal(size=t)
    loading = 1.0 + 0.5 * rng.normal(size=(n, 1))
    shocks = rng.normal(size=(n, t)) + factor * loading * common
    y = np.empty((n, t))
    y[:, 0] = shocks[:, 0]
    for k in range(1, t):
        y[:, k] = rho * y[:, k - 1] + shocks[:, k]
    y += rng.normal(scale=2.0, size=(n, 1))
    return pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n + 1), t),
            "t": np.tile(np.arange(1, t + 1), n),
            "y": y.ravel(),
        }
    )


def naive_cadf(frame: pd.DataFrame, lags: int, trend: bool) -> dict[int, float]:
    """Pesaran's unit t ratios written with pandas: the cross-section mean is a column."""
    data = frame.copy()
    data["m"] = data.groupby("t")["y"].transform("mean")
    out = {}
    for unit, g in data.groupby("id"):
        g = g.sort_values("t").reset_index(drop=True)
        cols = {
            "dy": g["y"].diff(),
            "ly": g["y"].shift(1),
            "lm": g["m"].shift(1),
            "dm": g["m"].diff(),
        }
        for j in range(1, lags + 1):
            cols[f"dy{j}"] = g["y"].diff().shift(j)
            cols[f"dm{j}"] = g["m"].diff().shift(j)
        if trend:
            cols["trend"] = g["t"].astype(float)
        design = pd.DataFrame(cols).dropna()
        x = design.drop(columns="dy")
        x.insert(0, "const", 1.0)
        order = ["const", "ly", *[c for c in x.columns if c not in ("const", "ly")]]
        matrix = x[order].to_numpy()
        response = design["dy"].to_numpy()
        beta, *_ = np.linalg.lstsq(matrix, response, rcond=None)
        resid = response - matrix @ beta
        sigma2 = resid @ resid / (len(response) - matrix.shape[1])
        out[int(unit)] = float(beta[1] / math.sqrt(sigma2 * np.linalg.inv(matrix.T @ matrix)[1, 1]))
    return out


@pytest.mark.parametrize(("name", "lags", "deterministics"), list(PLM))
def test_the_cips_statistic_is_the_one_of_plm_cipstest(
    name: str, lags: int, deterministics: str
) -> None:
    panel = declare(pd.read_csv(DATA / f"{name}.csv"))
    result = cips(panel, "y", lags=lags, deterministics=deterministics)  # type: ignore[arg-type]
    assert result.options["t_bar"] == pytest.approx(PLM[(name, lags, deterministics)], abs=1e-12)


@pytest.mark.parametrize(("lags", "trend"), [(0, False), (1, False), (2, True)])
def test_the_unit_statistics_follow_the_regression_written_directly(lags: int, trend: bool) -> None:
    frame = simulate(1, n=12, t=25)
    # drop the first periods of some units and the last of others: an unbalanced panel
    frame = frame[
        ~((frame["id"] <= 3) & (frame["t"] <= 2) | (frame["id"] == 5) & (frame["t"] > 22))
    ]
    result = cips(
        declare(frame.reset_index(drop=True)),
        "y",
        lags=lags,
        deterministics="trend" if trend else "constant",
        truncate=False,
    )
    assert result.details is not None
    expected = naive_cadf(frame, lags, trend)
    for unit, value in expected.items():
        assert result.details.loc[unit, "t_raw"] == pytest.approx(value, rel=1e-9)


def test_z_tbar_is_the_standardised_sum_with_the_tabulated_moments_and_a_left_tail_pvalue() -> None:
    result = cips(declare(simulate(2, n=30, t=30)), "y", lags=1)
    d = result.details
    assert d is not None
    z = (d["t"] - d["mean"]).sum() / math.sqrt((d["sd"] ** 2).sum())
    assert result.statistic == pytest.approx(z, rel=1e-12)
    assert result.pvalue == pytest.approx(stats.norm.cdf(z), rel=1e-12)
    assert result.statistic_name == "Z-t-bar" and result.n_groups == 30
    assert result.options["t_bar"] == pytest.approx(d["t"].mean(), rel=1e-12)
    # 30 units and 30 periods are tabulated as such: the moments are those of the (30, 30) cell
    assert (d["mean"] == d["mean"].iloc[0]).all() and d["mean"].iloc[0] == pytest.approx(-1.75)


def test_the_tables_are_read_at_the_first_key_at_or_above() -> None:
    # 40 units -> column 50, 17 periods -> row 20 (as pescadf does, without interpolation)
    result = cips(declare(simulate(3, n=40, t=17)), "y", lags=1)
    cell = result.details
    assert cell is not None
    assert cell["mean"].iloc[0] == pytest.approx(-1.73)  # T = 20, N = 50 of the constant case
    assert result.options["cv_5"] == pytest.approx(-2.11)  # truncated 5 % value at (20, 50)
    assert result.options["cv_1"] < result.options["cv_5"] < result.options["cv_10"]


def test_the_truncation_bounds_the_unit_statistics() -> None:
    panel = declare(simulate(4, n=30, t=30, rho=0.2))
    raw = cips(panel, "y", lags=1, truncate=False)
    cut = cips(panel, "y", lags=1)
    assert raw.details is not None and cut.details is not None
    assert (raw.details["t"] == raw.details["t_raw"]).all()
    assert cut.details["t"].min() >= -6.19 and cut.details["t"].max() <= 2.61
    np.testing.assert_allclose(cut.details["t"], raw.details["t_raw"].clip(-6.19, 2.61))
    trended = cips(panel, "y", lags=1, deterministics="trend")
    assert trended.details is not None
    assert trended.details["t"].min() >= -6.42 and trended.details["t"].max() <= 1.70


def test_a_constant_shift_changes_nothing_so_the_initial_cross_section_mean_has_no_option() -> None:
    # Pesaran advises removing the cross-section mean of the first period; it is a constant shift
    # of the series, which the constant of each regression absorbs (pescadf's nodemean gives the
    # same numbers, see the Stata cases).  A period effect added to the levels is not absorbed.
    frame = simulate(5)
    base = cips(declare(frame), "y", lags=1)
    shifted = frame.assign(y=frame["y"] - frame.loc[frame["t"] == 1, "y"].mean())
    assert cips(declare(shifted), "y", lags=1).statistic == pytest.approx(base.statistic, rel=1e-10)


def test_a_lag_for_each_unit_follows_the_order_of_the_units() -> None:
    panel = declare(simulate(6, n=10, t=30))
    per_unit = [0, 1, 2, 0, 1, 2, 0, 1, 2, 1]
    result = cips(panel, "y", lags=per_unit)
    assert result.details is not None
    assert list(result.details["lags"]) == [float(v) for v in per_unit]
    for unit, p in enumerate(per_unit, start=1):
        alone = cips(panel, "y", lags=p).details
        assert alone is not None
        assert result.details.loc[unit, "t_raw"] == pytest.approx(alone.loc[unit, "t_raw"])


def test_the_critical_values_need_a_balanced_panel() -> None:
    frame = simulate(7, n=20, t=30)
    balanced = cips(declare(frame), "y", lags=1)
    assert math.isfinite(balanced.options["cv_5"])  # type: ignore[arg-type]
    unbalanced = cips(declare(frame[~((frame["id"] == 4) & (frame["t"] > 25))]), "y", lags=1)
    assert math.isnan(unbalanced.options["cv_5"])  # type: ignore[arg-type]
    assert "balanced panels only" in unbalanced.summary("en")
    assert "balanced panels only" not in balanced.summary("en")
    assert "as if the unit statistics were independent" in balanced.summary("en")
    assert "comme si les statistiques des unités étaient" in balanced.summary("fr")


def test_it_rejects_stationary_series_and_not_series_with_a_common_factor_and_a_unit_root() -> None:
    assert cips(declare(simulate(8, n=40, t=40, rho=0.5)), "y", lags=1).rejects(0.01)
    assert not cips(declare(simulate(8, n=40, t=40, rho=1.0)), "y", lags=1).rejects(0.05)


def test_refusals() -> None:
    panel = declare(simulate(9, n=10, t=30))
    with pytest.raises(SpecificationError, match="deterministics must be one of"):
        cips(panel, "y", lags=1, deterministics="none")  # type: ignore[arg-type]
    for bad in (-1, 1.5, True, "1", []):
        with pytest.raises(SpecificationError, match="lags must be a non-negative integer"):
            cips(panel, "y", lags=bad)  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="gives 3 values for 10 units"):
        cips(panel, "y", lags=[1, 1, 1])
    with pytest.raises(NumericalError, match="no degree of freedom"):
        cips(declare(simulate(9, n=10, t=8)), "y", lags=1)
    twin = simulate(9, n=10, t=30)
    twin["y"] = np.tile(np.arange(30.0), 10)  # every unit equals the cross-section mean
    with pytest.raises(NumericalError, match="collinear"):
        cips(declare(twin), "y", lags=1)
    with pytest.raises(DataError):
        cips(declare(simulate(9, n=1, t=30)), "y", lags=1)
    with pytest.raises(DataError, match="up to 999 units"):
        cips(declare(simulate(9, n=1000, t=12)), "y", lags=0)
    gap = simulate(9, n=10, t=30)
    with pytest.raises(DataError):
        cips(declare(gap[~((gap["id"] == 2) & (gap["t"] == 10))]), "y", lags=1)


def test_by_name_alias_family_function_and_method_in_both_languages() -> None:
    panel = declare(simulate(10))
    direct = cips(panel, "y", lags=1)
    assert run_test("pescadf", panel, "y", lags=1).statistic == direct.statistic
    assert xt.unit_root(panel, "y", test="cips", lags=1).statistic == direct.statistic
    assert panel.unit_root("y", test="pesaran_cips", lags=1).pvalue == direct.pvalue
    assert "Pesaran CIPS panel unit-root test" in direct.summary("en")
    assert "Test de racine unitaire en panel CIPS de Pesaran" in direct.summary("fr")
