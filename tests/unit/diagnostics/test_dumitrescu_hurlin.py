"""Dumitrescu-Hurlin test: R's plm values, the finite-sample moments, a direct oracle, refusals."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import dumitrescu_hurlin, run_test
from pyxtpanel.diagnostics.causality.dumitrescu_hurlin import _unit_wald

DATA = Path(__file__).resolve().parents[2] / "golden" / "data"

#: plm 2.6.7, R 4.6.1: pgrangertest(y ~ x, data, index = c("id", "t"), order = K, test = ...) on the
#: frozen panels granger_yes and granger_no: (Wbar, Zbar, p of Zbar, Ztilde, p of Ztilde).
PLM = {
    ("granger_yes", 1): (
        5.894006816930,
        15.476208425890,
        5.02212782396e-54,
        13.173012705586,
        1.25488403253e-39,
    ),
    ("granger_yes", 2): (
        7.229524474739,
        11.693572215515,
        1.37479907473e-31,
        9.350565206761,
        8.71805874773e-21,
    ),
    ("granger_yes", 3): (
        8.375095444329,
        9.813536745341,
        9.84559732481e-23,
        7.231269979825,
        4.78498013427e-13,
    ),
    ("granger_no", 1): (
        1.526109095256,
        1.663703038741,
        0.0961717887273,
        1.212447853180,
        0.225340992516,
    ),
    ("granger_no", 2): (
        3.118357406466,
        2.500723183999,
        0.0123940011851,
        1.721796125935,
        0.0851064566356,
    ),
    ("granger_no", 3): (
        4.057349854778,
        1.930447888788,
        0.0535513661422,
        1.038438312073,
        0.299066040745,
    ),
}


def declare(frame: pd.DataFrame) -> xt.PanelData:
    return xt.xtset(frame, "id", "t")


def simulate(seed: int, *, n: int = 15, t: int = 30, slope: float = 0.0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    x = np.empty((n, t))
    y = np.empty((n, t))
    x[:, 0], y[:, 0] = rng.normal(size=n), rng.normal(size=n)
    for k in range(1, t):
        x[:, k] = 0.5 * x[:, k - 1] + rng.normal(size=n)
        y[:, k] = 0.4 * y[:, k - 1] + slope * x[:, k - 1] + rng.normal(size=n)
    return pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n + 1), t),
            "t": np.tile(np.arange(1, t + 1), n),
            "y": y.ravel(),
            "x": x.ravel(),
        }
    )


@pytest.mark.parametrize(("name", "lags"), list(PLM))
def test_the_statistics_are_those_of_plm_pgrangertest(name: str, lags: int) -> None:
    panel = declare(pd.read_csv(DATA / f"{name}.csv"))
    wbar, zbar, p_zbar, ztilde, p_ztilde = PLM[(name, lags)]
    assert dumitrescu_hurlin(
        panel, "y", "x", lags=lags, statistic="wbar"
    ).statistic == pytest.approx(wbar, rel=1e-11)
    z = dumitrescu_hurlin(panel, "y", "x", lags=lags, statistic="zbar")
    assert z.statistic == pytest.approx(zbar, rel=1e-11) and z.pvalue == pytest.approx(
        p_zbar, rel=1e-9
    )
    t = dumitrescu_hurlin(panel, "y", "x", lags=lags)
    assert t.statistic == pytest.approx(ztilde, rel=1e-11) and t.pvalue == pytest.approx(
        p_ztilde, rel=1e-9
    )


def test_the_unit_wald_statistic_follows_its_definition() -> None:
    frame = simulate(0)
    unit = frame[frame["id"] == 3]
    y, x = unit["y"].to_numpy(), unit["x"].to_numpy()
    k = 2
    response = y[k:]
    restricted = np.column_stack([np.ones(y.size - k), y[1:-1], y[:-2]])
    full = np.column_stack([restricted, x[1:-1], x[:-2]])

    def rss(design: np.ndarray) -> float:
        beta = np.linalg.solve(design.T @ design, design.T @ response)
        r = response - design @ beta
        return float(r @ r)

    expected = (rss(restricted) - rss(full)) / (rss(full) / (response.size - 2 * k - 1))
    assert _unit_wald(y, x, k, 3) == pytest.approx(expected, rel=1e-11)
    result = dumitrescu_hurlin(declare(frame), "y", "x", lags=k)
    assert result.details is not None
    assert result.details.loc[3, "wald"] == pytest.approx(expected, rel=1e-11)
    assert result.details.loc[3, "pvalue"] == pytest.approx(stats.chi2.sf(expected, k), rel=1e-12)
    assert result.n_obs == 15 * (30 - k) and result.n_groups == 15


def test_the_exact_finite_sample_moments_of_the_unit_statistic() -> None:
    # Under the null W_i has mean K (T-3K-1)/(T-3K-3) and variance 2K (T-3K-1)^2 (T-2K-3) /
    # ((T-3K-3)^2 (T-3K-5)); this is what Z-tilde standardises by.
    t, k, reps = 24, 2, 4000
    rng = np.random.default_rng(1)
    walds = np.empty(reps)
    for r in range(reps):
        y, x = rng.normal(size=t), rng.normal(size=t)
        walds[r] = _unit_wald(y, x, k, r)
    mean = k * (t - 3 * k - 1) / (t - 3 * k - 3)
    variance = (
        2 * k * (t - 3 * k - 1) ** 2 * (t - 2 * k - 3) / ((t - 3 * k - 3) ** 2 * (t - 3 * k - 5))
    )
    assert walds.mean() == pytest.approx(mean, rel=0.03)
    assert walds.var() == pytest.approx(variance, rel=0.12)


def test_it_rejects_when_x_causes_y_and_not_otherwise() -> None:
    assert dumitrescu_hurlin(declare(simulate(2, slope=0.6)), "y", "x").rejects(0.001)
    assert not dumitrescu_hurlin(declare(simulate(2)), "y", "x").rejects(0.001)


def test_the_choice_of_statistic_and_the_remark_for_wbar() -> None:
    panel = declare(simulate(3))
    wbar = dumitrescu_hurlin(panel, "y", "x", statistic="wbar")
    assert math.isnan(wbar.pvalue) and not wbar.rejects()
    assert "no p-value for W-bar" in wbar.summary("en")
    assert wbar.statistic_name == "W-bar" and wbar.options["wbar"] == wbar.statistic
    assert dumitrescu_hurlin(panel, "y", "x", statistic="zbar").statistic_name == "Z-bar"


def test_refusals() -> None:
    panel = declare(simulate(4))
    for bad in (0, -1, 1.5, True):
        with pytest.raises(SpecificationError, match="lags must be a positive integer"):
            dumitrescu_hurlin(panel, "y", "x", lags=bad)  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="statistic must be"):
        dumitrescu_hurlin(panel, "y", "x", statistic="w")  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="must be different variables"):
        dumitrescu_hurlin(panel, "y", "y")
    short = simulate(4, t=8)
    with pytest.raises(DataError):
        dumitrescu_hurlin(declare(short), "y", "x", lags=2)  # needs T > 3K + 5 = 11
    frame = simulate(4)
    frame["x"] = 1.0 + frame["y"].shift(1).fillna(0.0).where(frame["t"] > 1, 0.0)
    exact = frame.copy()
    exact["x"] = exact["y"]  # the lags of x are the lags of y: collinear
    with pytest.raises(NumericalError, match="collinear"):
        dumitrescu_hurlin(declare(exact), "y", "x")


def test_by_name_alias_family_function_and_method_in_both_languages() -> None:
    panel = declare(simulate(5))
    direct = dumitrescu_hurlin(panel, "y", "x", lags=2)
    assert run_test("granger", panel, "y", "x", lags=2).statistic == direct.statistic
    assert xt.causality(panel, "y", "x", test="xtgcause", lags=2).statistic == direct.statistic
    assert panel.causality("y", "x", test="dumitrescu_hurlin", lags=2).pvalue == direct.pvalue
    with pytest.raises(SpecificationError, match="choose a test"):
        xt.causality(panel, "y", "x")
    assert "Dumitrescu-Hurlin panel Granger non-causality test" in direct.summary("en")
    assert "Test de non-causalité de Granger en panel" in direct.summary("fr")
