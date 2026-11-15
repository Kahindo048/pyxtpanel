"""Transformations: within, between, two-way, quasi-demeaning, FOD, rolling (D.2)."""

import numpy as np
import pandas as pd
import pytest

from pyxtpanel import xtset
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.ops._spec import OPERATORS
from pyxtpanel.ops.within import twoway_demean

NAN = np.nan


def same(actual, expected, tol=1e-12):  # type: ignore[no-untyped-def]
    left, right = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    np.testing.assert_allclose(left, right, atol=tol, rtol=tol, equal_nan=True)


def unbalanced(seed: int = 0, n_units: int = 12, n_periods: int = 7) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for unit in range(1, n_units + 1):
        keep = rng.random(n_periods) < 0.75
        keep[rng.integers(n_periods)] = True
        for period in np.flatnonzero(keep):
            rows.append((unit, 2000 + int(period), rng.normal(), rng.normal()))
    frame = pd.DataFrame(rows, columns=["id", "t", "x", "z"])
    frame.loc[rng.random(len(frame)) < 0.1, "z"] = NAN
    return frame.sample(frac=1, random_state=seed).reset_index(drop=True)


@pytest.fixture
def panel():  # type: ignore[no-untyped-def]
    return xtset(unbalanced(), "id", "t")


# ---------------------------------------------------------------- within / between


def test_within_matches_pandas_groupby_and_xtdata(panel) -> None:  # type: ignore[no-untyped-def]
    frame = panel.to_frame()
    mean = frame.groupby("id")["x"].transform("mean")
    same(panel.within("x"), frame["x"] - mean)
    same(panel.within("x", add_mean=True), frame["x"] - mean + frame["x"].mean())
    same(panel.between("x"), mean)


def test_missing_values_use_the_variable_s_own_rows(panel) -> None:  # type: ignore[no-untyped-def]
    frame = panel.to_frame()
    expect = frame["z"] - frame.groupby("id")["z"].transform("mean")
    same(panel.within("z"), expect)
    assert panel.within("z").isna().equals(frame["z"].isna())


def test_listwise_drops_a_row_for_every_variable(panel) -> None:  # type: ignore[no-untyped-def]
    out = panel.within(["x", "z"], listwise=True)
    assert out["x"].isna().equals(panel.to_frame()["z"].isna())  # x lost the rows z lacks
    frame = panel.to_frame()
    keep = frame["z"].notna()
    expect = frame.loc[keep, "x"] - frame.loc[keep].groupby("id")["x"].transform("mean")
    same(out.loc[keep, "x"], expect)


def test_within_of_a_constant_is_zero_and_between_is_the_constant() -> None:
    frame = pd.DataFrame({"id": [1, 1, 2, 2], "t": [1, 2, 1, 2], "k": [3.0, 3, 5, 5]})
    panel = xtset(frame, "id", "t")
    same(panel.within("k"), [0, 0, 0, 0])
    same(panel.between("k"), [3, 3, 5, 5])


def test_group_statistics_by_unit_and_period(panel) -> None:  # type: ignore[no-untyped-def]
    frame = panel.to_frame()
    same(panel.group_mean("x", "time"), frame.groupby("t")["x"].transform("mean"))
    same(
        panel.group_count("z"),
        frame.groupby("id")["z"].transform("count").where(frame["z"].notna()),
    )
    same(panel.group_size("unit"), frame.groupby("id")["x"].transform("size"))
    same(panel.group_size("time"), frame.groupby("t")["x"].transform("size"))
    assert panel.group_mean("z").isna().equals(frame["z"].isna())
    assert not panel.group_mean("z", at_missing=True).isna().any()  # every unit has some z


def test_a_frame_of_results_keeps_the_names_and_the_index() -> None:
    frame = pd.DataFrame(
        {"id": [1, 1], "t": [1, 2], "x": [1.0, 3], "y": [2.0, 6]}, index=["a", "b"]
    )
    out = xtset(frame, "id", "t").within(["x", "y"])
    assert list(out.columns) == ["x", "y"] and out.index.tolist() == ["a", "b"]
    same(out["y"], [-2, 2])


# ---------------------------------------------------------------- two-way


def dummy_residual(frame: pd.DataFrame, column: str) -> np.ndarray:
    """Residual of the regression of ``column`` on unit and period dummies (reference)."""
    use = frame[column].notna().to_numpy()
    sub = frame.loc[use]
    design = pd.get_dummies(sub["id"], dtype=float)
    time = pd.get_dummies(sub["t"], dtype=float).iloc[:, 1:]
    matrix = np.column_stack([design.to_numpy(), time.to_numpy()])
    beta, *_ = np.linalg.lstsq(matrix, sub[column].to_numpy(), rcond=None)
    out = np.full(len(frame), NAN)
    out[use] = sub[column].to_numpy() - matrix @ beta
    return out


def test_twoway_unbalanced_equals_the_dummy_regression_residual(panel) -> None:  # type: ignore[no-untyped-def]
    frame = panel.to_frame()
    result = panel.within(["x", "z"], effects="twoway")
    same(result["x"], dummy_residual(frame, "x"), 1e-9)
    same(result["z"], dummy_residual(frame, "z"), 1e-9)
    report = result.attrs["convergence"]
    assert report["converged"] == (True, True) and report["exact"] == (False, False)


def test_twoway_closed_form_in_a_balanced_panel() -> None:
    rng = np.random.default_rng(3)
    frame = pd.DataFrame(
        {"id": np.repeat(np.arange(6), 5), "t": np.tile(np.arange(5), 6), "x": rng.normal(size=30)}
    )
    panel = xtset(frame, "id", "t")
    result = panel.within("x", effects="twoway")
    expect = (
        frame["x"]
        - frame.groupby("id")["x"].transform("mean")
        - frame.groupby("t")["x"].transform("mean")
        + frame["x"].mean()
    )
    same(result, expect)
    assert result.attrs["convergence"]["exact"] == (True,)
    assert result.attrs["convergence"]["iterations"] == (0,)
    same(result, dummy_residual(frame, "x"), 1e-10)


def test_twoway_effects_time_and_unit_only(panel) -> None:  # type: ignore[no-untyped-def]
    frame = panel.to_frame()
    same(panel.within("x", effects="time"), frame["x"] - frame.groupby("t")["x"].transform("mean"))


def test_twoway_does_not_converge_with_too_few_iterations(panel) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(NumericalError, match="did not converge") as info:
        panel.within("x", effects="twoway", maxiter=2, tol=1e-14)
    assert info.value.key == "ops.no_convergence"
    values = np.arange(10.0)
    _, report = twoway_demean(panel.index, np.resize(values, panel.n_obs), maxiter=2, strict=False)
    assert report.converged == (False,)


def test_twoway_on_a_column_with_a_dropped_unit_stays_exact() -> None:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 2, 2, 3, 3],
            "t": [1, 2, 1, 2, 1, 2],
            "x": [1.0, 2, 4, 3, NAN, NAN],  # unit 3 is entirely missing
        }
    )
    result = xtset(frame, "id", "t").within("x", effects="twoway")
    assert result.isna().tolist() == [False, False, False, False, True, True]
    assert result.attrs["convergence"]["exact"] == (True,)


# ---------------------------------------------------------------- quasi-demeaning


def test_quasi_demeaning_with_scalar_and_per_unit_theta(panel) -> None:  # type: ignore[no-untyped-def]
    frame = panel.to_frame()
    mean = frame.groupby("id")["x"].transform("mean")
    same(panel.quasi_demean("x", 1.0), frame["x"] - mean)
    same(panel.quasi_demean("x", 0.0), frame["x"])
    theta = panel.theta_re(sigma_u=1.3, sigma_e=0.7)
    assert theta.index.tolist() == panel.unit_labels.tolist() and ((theta > 0) & (theta < 1)).all()
    per_row = frame["id"].map(theta)
    same(panel.quasi_demean("x", theta.to_numpy()), frame["x"] - per_row * mean)


def test_theta_formula_and_validation(panel) -> None:  # type: ignore[no-untyped-def]
    t = panel.counts
    same(panel.theta_re(2.0, 1.0).to_numpy(), 1 - 1 / np.sqrt(1 + t * 4.0))
    with pytest.raises(SpecificationError, match="theta"):
        panel.quasi_demean("x", [0.1, 0.2])
    with pytest.raises(SpecificationError):
        panel.theta_re(-1.0, 1.0)


# ---------------------------------------------------------------- forward orthogonal deviations


def fod_reference(values: list[float], scale: bool = True) -> list[float]:
    out = []
    for i, v in enumerate(values):
        later = [w for w in values[i + 1 :] if not np.isnan(w)]
        if np.isnan(v) or not later:
            out.append(NAN)
        else:
            f = len(later)
            out.append((np.sqrt(f / (f + 1)) if scale else 1.0) * (v - np.mean(later)))
    return out


def test_fod_follows_the_next_observed_rows_not_the_dates() -> None:
    frame = pd.DataFrame(
        {"id": [1] * 5 + [2] * 2, "t": [1, 2, 4, 5, 9, 1, 2], "x": [1.0, 2, NAN, 7, 11, 5, 6]}
    )
    panel = xtset(frame, "id", "t")
    got = panel.fod("x")
    same(got.iloc[:5], fod_reference([1.0, 2, NAN, 7, 11]))
    same(got.iloc[5:], fod_reference([5.0, 6]))
    same(panel.fod("x", scale=False).iloc[:5], fod_reference([1.0, 2, NAN, 7, 11], scale=False))
    assert np.isnan(got.iloc[4]) and np.isnan(got.iloc[2]) and np.isnan(got.iloc[6])


def test_fod_on_random_unbalanced_panels(panel) -> None:  # type: ignore[no-untyped-def]
    sorted_frame = panel.to_frame(order="sorted")
    expect = []
    for _, group in sorted_frame.groupby("id", sort=True):
        expect += fod_reference(group["z"].tolist())
    out = panel.fod("z")
    same(out.loc[sorted_frame.index], expect)


def test_fod_of_a_constant_plus_noise_is_translation_invariant(panel) -> None:  # type: ignore[no-untyped-def]
    shifted = panel.to_frame().assign(big=lambda d: d["x"] + 1e9)
    other = xtset(shifted, "id", "t")
    same(other.fod("big"), panel.fod("x"), 1e-6)


# ---------------------------------------------------------------- rolling, cumulative, growth


def test_rolling_respects_calendar_gaps() -> None:
    frame = pd.DataFrame({"id": [1] * 5, "t": [1, 2, 4, 5, 6], "x": [1.0, 2, 4, 5, 6]})
    panel = xtset(frame, "id", "t")
    same(panel.rolling("x", 2, "mean"), [NAN, 1.5, NAN, 4.5, 5.5])  # period 3 is not in any window
    same(panel.rolling("x", 3, "sum", min_obs=2), [NAN, 3.0, 6.0, 9.0, 15.0])
    same(panel.rolling("x", 2, "std"), [NAN, np.sqrt(0.5), NAN, np.sqrt(0.5), np.sqrt(0.5)])
    same(panel.rolling("x", 3, "max", min_obs=1), [1.0, 2.0, 4.0, 5.0, 6.0])
    same(panel.rolling("x", 2, "min", min_obs=1), [1.0, 1.0, 4.0, 4.0, 5.0])


def test_rolling_validation() -> None:
    panel = xtset(pd.DataFrame({"id": [1, 1], "t": [1, 2], "x": [1.0, 2]}), "id", "t")
    for window, min_obs in ((0, None), (2, 3), (2, 0)):
        with pytest.raises(SpecificationError, match="window"):
            panel.rolling("x", window, min_obs=min_obs)
    with pytest.raises(SpecificationError, match="unknown statistic"):
        panel.rolling("x", 2, "median")  # type: ignore[arg-type]


def test_cumulative_and_growth() -> None:
    frame = pd.DataFrame({"id": [1, 1, 1, 2, 2], "t": [1, 2, 3, 1, 2], "x": [1.0, NAN, 3, 4, 6]})
    panel = xtset(frame, "id", "t")
    same(panel.cumulative("x"), [1.0, NAN, 4.0, 4.0, 10.0])
    same(panel.cumulative("x", "mean"), [1.0, NAN, 2.0, 4.0, 5.0])
    same(panel.growth("x"), [NAN, NAN, NAN, NAN, 0.5])
    same(panel.growth("x", log=True), [NAN, NAN, NAN, NAN, np.log(1.5)])
    zero = xtset(pd.DataFrame({"id": [1, 1], "t": [1, 2], "x": [0.0, 1]}), "id", "t")
    assert zero.growth("x").isna().all()  # division by zero is missing, not infinite


# ---------------------------------------------------------------- generic access and registry


def test_transform_by_name_and_registry_provenance(panel) -> None:  # type: ignore[no-untyped-def]
    same(panel.transform("within", "x"), panel.within("x"))
    same(panel.transform("rolling", "x", window=2, stat="sum"), panel.rolling("x", 2, "sum"))
    expected = {"between", "within", "quasi_demean", "twoway", "fod", "rolling", "cumulative"}
    assert expected <= set(OPERATORS.names())
    assert all(OPERATORS.spec(n).summary for n in OPERATORS.names())
    with pytest.raises(SpecificationError, match="unknown operator"):
        panel.transform("nope", "x")
