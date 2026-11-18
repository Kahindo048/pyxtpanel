"""Westerlund test: a direct implementation, the moments table, invariances, refusals."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import pyxtpanel as xt
from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import run_test, westerlund
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


def by_hand(
    frame: pd.DataFrame, names: list[str], trend: bool, alternative: str, mean: float, var: float
) -> float:
    """The statistic written with loops over units and periods."""
    nums, sqs, lengths = [], [], []
    for _, g in frame.groupby("id"):
        t = len(g)
        design = [np.ones(t)] + ([np.arange(1.0, t + 1)] if trend else [])
        design += [g[c].to_numpy() for c in names]
        x = np.column_stack(design)
        y = g["y"].to_numpy()
        e = y - x @ np.linalg.solve(x.T @ x, x.T @ y)
        total = 0.0
        cumulated = 0.0
        for value in e:
            cumulated += value
            total += cumulated**2
        nums.append(total)
        sqs.append(float(e @ e))
        lengths.append(t)
    ratio = (
        float(np.mean([a / b for a, b in zip(nums, sqs, strict=True)]))
        if alternative == "some"
        else sum(nums) / sum(sqs)
    )
    return math.sqrt(len(nums)) * (ratio / max(lengths) ** 2 - mean) / math.sqrt(var)


@pytest.mark.parametrize("alternative", ["some", "all"])
@pytest.mark.parametrize("trend", [False, True])
def test_the_statistic_follows_its_definition_with_the_tabulated_moments(
    alternative: str, trend: bool
) -> None:
    frame = simulate(0, k=2)
    # drop the first periods of two units and the last of another: an unbalanced panel
    early = frame["id"].isin([2, 3]) & (frame["t"] <= 4)
    late = (frame["id"] == 5) & (frame["t"] > 35)
    frame = frame[~(early | late)]
    table = load_table("westerlund_moments.csv")
    row = table[
        (table["m"] == 2)
        & (table["deterministics"] == ("trend" if trend else "constant"))
        & (table["alternative"] == alternative)
    ].iloc[0]
    expected = by_hand(frame, ["x1", "x2"], trend, alternative, row["mean"], row["variance"])
    result = westerlund(
        declare(frame),
        "y",
        ["x1", "x2"],
        alternative=alternative,  # type: ignore[arg-type]
        deterministics="trend" if trend else "constant",
    )
    assert result.statistic == pytest.approx(expected, rel=1e-11)
    assert result.pvalue == pytest.approx(stats.norm.cdf(expected), rel=1e-11)
    assert result.statistic_name == "VR" and result.n_groups == 20


def test_the_moments_table_has_every_cell_and_decreasing_means() -> None:
    table = load_table("westerlund_moments.csv")
    assert len(table) == 7 * 2 * 2 and (table["variance"] > 0).all()
    for _, group in table.groupby(["deterministics", "alternative"]):
        assert group.sort_values("m")["mean"].is_monotonic_decreasing
    # a trend lowers the mean for the same number of regressors
    wide = table.pivot_table(index=["m", "alternative"], columns="deterministics", values="mean")
    assert (wide["trend"] < wide["constant"]).all()


def test_it_rejects_cointegrated_panels_and_not_spurious_ones() -> None:
    cointegrated = declare(simulate(1, cointegrated=True))
    assert westerlund(cointegrated, "y", ["x1"]).rejects(0.001)
    assert westerlund(cointegrated, "y", ["x1"], alternative="all").rejects(0.001)
    assert not westerlund(declare(simulate(1)), "y", ["x1"]).rejects(0.001)


def test_unit_effects_and_the_scale_do_not_change_the_statistic() -> None:
    frame = simulate(2)
    base = westerlund(declare(frame), "y", ["x1"])
    shifted = frame.assign(y=frame["y"] + 5.0 * frame["id"])  # a different constant per unit
    scaled = frame.assign(y=1000 * frame["y"], x1=1000 * frame["x1"])
    for other in (shifted, scaled):
        assert westerlund(declare(other), "y", ["x1"]).statistic == pytest.approx(
            base.statistic, rel=1e-7
        )


def test_a_trend_of_each_unit_is_absorbed_with_deterministics_trend_only() -> None:
    frame = simulate(2)
    trended = frame.assign(y=frame["y"] + 0.3 * frame["t"] * frame["id"])
    plain = westerlund(declare(frame), "y", ["x1"], deterministics="trend")
    absorbed = westerlund(declare(trended), "y", ["x1"], deterministics="trend")
    assert absorbed.statistic == pytest.approx(plain.statistic, rel=1e-6)
    assert westerlund(declare(trended), "y", ["x1"]).statistic != pytest.approx(
        westerlund(declare(frame), "y", ["x1"]).statistic
    )


def test_demean_removes_the_period_means_of_every_variable() -> None:
    frame = simulate(3, k=2)
    by_period = frame.copy()
    for column in ("y", "x1", "x2"):
        by_period[column] -= by_period.groupby("t")[column].transform("mean")
    direct = westerlund(declare(by_period), "y", ["x1", "x2"])
    assert westerlund(declare(frame), "y", ["x1", "x2"], demean=True).statistic == pytest.approx(
        direct.statistic, rel=1e-10
    )


def test_t_is_the_longest_unit_on_an_unbalanced_panel() -> None:
    frame = simulate(4)
    cut = frame[~((frame["id"] <= 10) & (frame["t"] > 25))]
    result = westerlund(declare(cut), "y", ["x1"])
    assert result.details is not None
    assert result.details["periods"].max() == 40 and result.details["periods"].min() == 25
    # every unit ratio is divided by the square of the longest unit (40), whatever its length
    assert result.options["variance_ratio"] == pytest.approx(
        float(result.details["variance_ratio"].mean()), rel=1e-12
    )


def test_refusals() -> None:
    panel = declare(simulate(5, k=2))
    for bad in ([], ["y"], ["x1", "x1"], [f"x{i}" for i in range(8)]):
        with pytest.raises(SpecificationError, match="between 1 and 7"):
            westerlund(panel, "y", bad)
    for name, value in (("alternative", "most"), ("deterministics", "none")):
        with pytest.raises(SpecificationError, match=f"{name} must be one of"):
            westerlund(panel, "y", ["x1"], **{name: value})  # type: ignore[arg-type]
    frame = simulate(5)
    with pytest.raises(DataError):
        westerlund(declare(frame[~((frame["id"] == 3) & (frame["t"] == 10))]), "y", ["x1"])
    with pytest.raises(DataError):
        westerlund(declare(frame[frame["t"] <= 3]), "y", ["x1"], deterministics="trend")
    exact = frame.assign(y=2.0 * frame["x1"] + 3.0)
    with pytest.raises(NumericalError, match="fits perfectly"):
        westerlund(declare(exact), "y", ["x1"])


def test_by_name_family_function_and_method_in_both_languages() -> None:
    panel = declare(simulate(6))
    direct = westerlund(panel, "y", ["x1"])
    assert run_test("westerlund", panel, "y", ["x1"]).statistic == direct.statistic
    assert xt.cointegration(panel, "y", ["x1"], test="westerlund").statistic == direct.statistic
    assert panel.cointegration("y", "x1", test="westerlund").pvalue == direct.pvalue
    assert "Westerlund variance-ratio test for cointegration" in direct.summary("en")
    assert "Test de cointégration de Westerlund" in direct.summary("fr")
