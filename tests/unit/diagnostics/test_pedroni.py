"""Pedroni tests: a direct implementation, the adjustment table, invariances, refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import pedroni, run_test
from pyxtpanel.diagnostics._shared.tables import load_table


def simulate(
    seed: int, *, n: int = 20, t: int = 40, cointegrated: bool = False, k: int = 1
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    x = np.cumsum(rng.normal(size=(n, t, k)), axis=1)
    e = rng.normal(size=(n, t)) if cointegrated else np.cumsum(rng.normal(size=(n, t)), axis=1)
    y = rng.normal(scale=2.0, size=(n, 1)) + x.sum(axis=2) + e
    frame = {
        "id": np.repeat(np.arange(1, n + 1), t),
        "t": np.tile(np.arange(1, t + 1), n),
        "y": y.ravel(),
    }
    for j in range(k):
        frame[f"x{j + 1}"] = x[:, :, j].ravel()
    return pd.DataFrame(frame)


def declare(frame: pd.DataFrame) -> xt.PanelData:
    return xt.xtset(frame.reset_index(drop=True), "id", "t")


def long_run(w: np.ndarray, lags: int) -> np.ndarray:
    """Bartlett long-run covariance of the centred series (weights 1 - j / (lags + 1))."""
    w = w - w.mean(axis=0)
    n = len(w)
    omega = w.T @ w / n
    for j in range(1, lags + 1):
        gamma = w[j:].T @ w[:-j] / n
        omega = omega + (1 - j / (lags + 1)) * (gamma + gamma.T)
    return omega


def by_hand(frame: pd.DataFrame, names: list[str], lags: int, p: int) -> dict[str, float]:
    """The seven raw statistics written unit by unit (constant case, Bartlett kernel)."""
    parts = []
    for _, g in frame.groupby("id"):
        y = g["y"].to_numpy()
        x = g[names].to_numpy()
        t = len(y)
        design = np.column_stack([np.ones(t), x])
        e = y - design @ np.linalg.solve(design.T @ design, design.T @ y)
        w = np.column_stack([np.diff(y), np.diff(x, axis=0)])
        om = long_run(w, lags)
        l11 = om[0, 0] - om[0, 1:] @ np.linalg.solve(om[1:, 1:], om[1:, 0])
        lag, cur = e[:-1], e[1:]
        rho = (lag @ cur) / (lag @ lag)
        u = cur - rho * lag
        s2 = float(np.mean((u - u.mean()) ** 2))
        sig2 = float(long_run(u[:, None], lags)[0, 0])
        num = float(lag @ (cur - lag)) - len(lag) * 0.5 * (sig2 - s2)
        de = np.diff(e)
        rows = range(p + 1, t)
        yy = np.array([de[r - 1] for r in rows])
        xx = np.array([[e[r - 1], *[de[r - 1 - j] for j in range(1, p + 1)]] for r in rows])
        ua = yy - xx @ np.linalg.lstsq(xx, yy, rcond=None)[0]
        s2a = float(np.mean((ua - ua.mean()) ** 2))
        parts.append(
            (
                l11,
                float(lag @ lag),
                num,
                sig2,
                float(xx[:, 0] @ yy),
                float(xx[:, 0] @ xx[:, 0]),
                s2a,
            )
        )
    l11, sxx, num, sig2, na, sa, s2a = (np.array(c) for c in zip(*parts, strict=True))
    n, t = len(parts), frame.groupby("id").size().max() - 1.0
    ws, wn = np.sum(sxx / l11), np.sum(num / l11)
    return {
        "panel_nu": t * t * n**1.5 / ws,
        "panel_rho": t * math.sqrt(n) * wn / ws,
        "panel_pp": wn / math.sqrt(np.mean(sig2 / l11) * ws),
        "panel_adf": np.sum(na / l11) / math.sqrt(np.mean(s2a / l11) * np.sum(sa / l11)),
        "group_rho": t / math.sqrt(n) * np.sum(num / sxx),
        "group_pp": np.sum(num / np.sqrt(sig2 * sxx)) / math.sqrt(n),
        "group_adf": np.sum(na / np.sqrt(s2a * sa)) / math.sqrt(n),
    }


def standardised(raw: float, label: str, n: int, m: int) -> float:
    table = load_table("pedroni_moments.csv")
    row = table[
        (table["deterministics"] == "constant") & (table["m"] == m) & (table["statistic"] == label)
    ].iloc[0]
    return (raw - row["mean"] * math.sqrt(n)) / math.sqrt(row["variance"])


@pytest.mark.parametrize(("lags", "p"), [(0, 1), (2, 1), (3, 2)])
def test_the_seven_statistics_follow_their_definitions_on_an_unbalanced_panel(
    lags: int, p: int
) -> None:
    frame = simulate(0, k=2)
    early = frame["id"].isin([2, 3]) & (frame["t"] <= 5)
    late = (frame["id"] == 6) & (frame["t"] > 33)
    frame = frame[~(early | late)]
    raw = by_hand(frame, ["x1", "x2"], lags, p)
    panel = declare(frame)
    for ar, keys in (("same", ("nu", "rho", "pp", "adf")), ("panelspecific", ("rho", "pp", "adf"))):
        for key in keys:
            result = pedroni(
                panel,
                "y",
                ["x1", "x2"],
                kernel_lags=lags,
                adf_lags=p,
                ar=ar,  # type: ignore[arg-type]
                statistic=key,  # type: ignore[arg-type]
            )
            label = f"{'panel' if ar == 'same' else 'group'}_{key}"
            expected = standardised(raw[label], label, 20, 2)
            assert result.statistic == pytest.approx(expected, rel=1e-10, abs=1e-10)
            tail = stats.norm.sf(expected) if key == "nu" else stats.norm.cdf(expected)
            assert result.pvalue == pytest.approx(tail, rel=1e-9)


def test_the_adjustment_table_has_every_cell_and_the_known_corners() -> None:
    table = load_table("pedroni_moments.csv")
    assert len(table) == 3 * 7 * 7 and (table["variance"] > 0).all()
    corner = table[
        (table["deterministics"] == "constant")
        & (table["m"] == 1)
        & (table["statistic"] == "panel_nu")
    ].iloc[0]
    assert (corner["mean"], corner["variance"]) == (8.62, 60.75)
    # the printed Table 2 of Pedroni (1999) for two regressors with heterogeneous intercepts
    row = table[
        (table["deterministics"] == "constant")
        & (table["m"] == 2)
        & (table["statistic"] == "group_pp")
    ].iloc[0]
    assert (row["mean"], row["variance"]) == (-2.453, 0.618)
    assert (table[table["statistic"] == "panel_nu"]["mean"] > 0).all()
    others = table[table["statistic"] != "panel_nu"]
    positive = others[others["mean"] >= 0]
    # Stata's table has +25.889 where the series and Pedroni's table have -25.889 (a sign error)
    assert positive[["deterministics", "m", "statistic"]].values.tolist() == [
        ["none", 7, "panel_rho"]
    ]


def test_it_rejects_cointegrated_panels_and_not_spurious_ones() -> None:
    cointegrated = declare(simulate(1, cointegrated=True))
    spurious = declare(simulate(1))
    for ar, key in (
        ("same", "adf"),
        ("same", "nu"),
        ("panelspecific", "adf"),
        ("panelspecific", "pp"),
    ):
        kwargs = {"kernel_lags": 2, "ar": ar, "statistic": key}
        assert pedroni(cointegrated, "y", ["x1"], **kwargs).rejects(0.001)  # type: ignore[arg-type]
        assert not pedroni(spurious, "y", ["x1"], **kwargs).rejects(0.001)  # type: ignore[arg-type]


def test_details_list_the_family_with_one_sided_p_values() -> None:
    panel = declare(simulate(2))
    same = pedroni(panel, "y", ["x1"], kernel_lags=1, ar="same")
    group = pedroni(panel, "y", ["x1"], kernel_lags=1)
    assert same.details is not None and group.details is not None
    assert list(same.details.index) == [
        "Panel variance ratio",
        "Panel modified Phillips-Perron t",
        "Panel Phillips-Perron t",
        "Panel ADF t",
    ]
    assert len(group.details) == 3 and group.statistic_name == "Group ADF t"
    z = same.details["statistic"].to_numpy()
    p = same.details["pvalue"].to_numpy()
    assert p[0] == pytest.approx(stats.norm.sf(z[0])) and p[1:] == pytest.approx(
        stats.norm.cdf(z[1:])
    )
    assert same.options["adf_lags"] == 1 and same.n_groups == 20 and same.n_obs == 20 * 39


def test_the_sign_error_of_the_stata_table_is_flagged() -> None:
    panel = declare(simulate(9, k=7))
    flagged = pedroni(
        panel, "y", [f"x{j}" for j in range(1, 8)], kernel_lags=1, ar="same", deterministics="none"
    )
    assert "sign looks wrong in Stata" in flagged.summary("en")
    plain = pedroni(panel, "y", ["x1", "x2"], kernel_lags=1, ar="same", deterministics="none")
    assert "sign looks wrong" not in plain.summary("en")


def test_unit_effects_scale_and_trends_follow_the_deterministic_terms() -> None:
    frame = simulate(3)
    base = pedroni(declare(frame), "y", ["x1"], kernel_lags=2, ar="same")
    shifted = frame.assign(y=frame["y"] + 5.0 * frame["id"])  # constants are absorbed
    scaled = frame.assign(y=1000 * frame["y"], x1=1000 * frame["x1"])
    for other in (shifted, scaled):
        assert pedroni(
            declare(other), "y", ["x1"], kernel_lags=2, ar="same"
        ).statistic == pytest.approx(base.statistic, rel=1e-6)
    trended = frame.assign(y=frame["y"] + 0.3 * frame["t"] * frame["id"])
    kwargs = {"kernel_lags": 2, "ar": "same", "deterministics": "trend"}
    assert pedroni(declare(trended), "y", ["x1"], **kwargs).statistic == pytest.approx(  # type: ignore[arg-type]
        pedroni(declare(frame), "y", ["x1"], **kwargs).statistic,
        rel=1e-6,  # type: ignore[arg-type]
    )


def test_demean_removes_the_period_means_of_every_variable() -> None:
    frame = simulate(4, k=2)
    by_period = frame.copy()
    for column in ("y", "x1", "x2"):
        by_period[column] -= by_period.groupby("t")[column].transform("mean")
    direct = pedroni(declare(by_period), "y", ["x1", "x2"], kernel_lags=2)
    assert pedroni(
        declare(frame), "y", ["x1", "x2"], kernel_lags=2, demean=True
    ).statistic == pytest.approx(direct.statistic, rel=1e-10)


def test_the_choice_of_augmented_lags_picks_a_number_in_range() -> None:
    panel = declare(simulate(5))
    for criterion in ("aic", "bic", "hqic"):
        result = pedroni(
            panel,
            "y",
            ["x1"],
            kernel_lags=2,
            adf_lags=4,
            adf_lag_selection=criterion,  # type: ignore[arg-type]
        )
        assert 1 <= int(result.options["adf_lags"]) <= 4  # type: ignore[call-overload]


def test_refusals() -> None:
    panel = declare(simulate(6, k=2))
    for bad in ([], ["y"], ["x1", "x1"], [f"x{i}" for i in range(8)]):
        with pytest.raises(SpecificationError, match="between 1 and 7"):
            pedroni(panel, "y", bad, kernel_lags=1)
    for value in (-1, 1.5, True, "2"):
        with pytest.raises(SpecificationError, match="kernel_lags must be a non-negative integer"):
            pedroni(panel, "y", ["x1"], kernel_lags=value)  # type: ignore[arg-type]
    for value in (0, -1, 1.5, True):
        with pytest.raises(SpecificationError, match="adf_lags must be an integer of at least 1"):
            pedroni(panel, "y", ["x1"], kernel_lags=1, adf_lags=value)  # type: ignore[arg-type]
    for name, value in (
        ("ar", "pooled"),
        ("deterministics", "drift"),
        ("kernel", "tukey"),
        ("adf_lag_selection", "ic"),
    ):
        with pytest.raises(SpecificationError, match=f"{name} must be one of"):
            pedroni(panel, "y", ["x1"], kernel_lags=1, **{name: value})  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="the statistic must be one of 'rho', 'pp', 'adf'"):
        pedroni(panel, "y", ["x1"], kernel_lags=1, statistic="nu")
    frame = simulate(6)
    with pytest.raises(DataError):
        pedroni(
            declare(frame[~((frame["id"] == 3) & (frame["t"] == 10))]), "y", ["x1"], kernel_lags=1
        )
    with pytest.raises(DataError):
        pedroni(declare(frame[frame["t"] <= 4]), "y", ["x1"], kernel_lags=1)
    exact = frame.assign(y=2.0 * frame["x1"] + 3.0)
    with pytest.raises(NumericalError, match="degenerate"):
        pedroni(declare(exact), "y", ["x1"], kernel_lags=1)


def test_the_kernel_lags_are_required_and_the_automatic_bandwidth_is_not_offered() -> None:
    panel = declare(simulate(7))
    with pytest.raises(SpecificationError, match="pedroni: missing a required keyword-only"):
        xt.cointegration(panel, "y", ["x1"], test="pedroni")
    with pytest.raises(SpecificationError, match="automatic bandwidth"):
        pedroni(panel, "y", ["x1"], kernel_lags=None)  # type: ignore[arg-type]


def test_by_name_family_function_and_method_in_both_languages() -> None:
    panel = declare(simulate(8))
    direct = pedroni(panel, "y", ["x1"], kernel_lags=2)
    assert run_test("pedroni", panel, "y", ["x1"], kernel_lags=2).statistic == direct.statistic
    assert (
        xt.cointegration(panel, "y", ["x1"], test="pedroni", kernel_lags=2).statistic
        == direct.statistic
    )
    assert panel.cointegration("y", "x1", test="pedroni", kernel_lags=2).pvalue == direct.pvalue
    assert "Pedroni panel cointegration test" in direct.summary("en")
    assert "Test de cointégration en panel de Pedroni" in direct.summary("fr")
