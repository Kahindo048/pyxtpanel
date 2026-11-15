"""RESET test: a nested F computed by hand as oracle, collinearity handling and refusals."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import reset, run_test


def data(seed: int = 0, n: int = 400, curvature: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 2))
    y = 1 + x @ np.array([1.0, -0.5]) + curvature * x[:, 0] ** 2 + rng.normal(size=n)
    return y, x


def nested_f(y: np.ndarray, base: np.ndarray, extra: np.ndarray) -> tuple[float, int, int]:
    """F test of the extra columns by two least-squares fits with numpy."""

    def rss(z: np.ndarray) -> tuple[float, int]:
        beta, _, rank, _ = np.linalg.lstsq(z, y, rcond=None)
        return float(np.sum((y - z @ beta) ** 2)), int(rank)

    rss_r, rank_r = rss(base)
    rss_u, rank_u = rss(np.column_stack([base, extra]))
    q, dfr = rank_u - rank_r, y.size - rank_u
    return ((rss_r - rss_u) / q) / (rss_u / dfr), q, dfr


def test_the_fitted_values_version_is_the_nested_f_of_the_raw_powers() -> None:
    y, x = data(curvature=0.4)
    z = np.column_stack([np.ones(y.size), x])
    fitted = z @ np.linalg.lstsq(z, y, rcond=None)[0]
    expected, q, dfr = nested_f(y, z, np.column_stack([fitted**2, fitted**3, fitted**4]))
    result = reset(y, x)
    assert result.statistic == pytest.approx(expected, rel=1e-8)
    assert result.df == (q, dfr) == (3.0, 394.0)
    assert result.pvalue == pytest.approx(float(stats.f.sf(expected, q, dfr)), rel=1e-7)


def test_the_rhs_version_adds_the_powers_of_every_regressor() -> None:
    y, x = data(curvature=0.4)
    z = np.column_stack([np.ones(y.size), x])
    powers = np.column_stack([x**2, x**3, x**4])
    expected, q, dfr = nested_f(y, z, powers)
    result = reset(y, x, rhs=True)
    assert result.statistic == pytest.approx(expected, rel=1e-8)
    assert result.df == (q, dfr) == (6.0, 391.0)
    assert result.options == {"rhs": True}


def test_a_binary_regressor_is_not_raised_to_powers_and_a_constant_one_is_dropped() -> None:
    y, x = data(curvature=0.4)
    dummy = (x[:, 1] > 0).astype(float)
    result = reset(y, np.column_stack([x[:, 0], dummy, np.ones(y.size)]), rhs=True)
    assert result.df[0] == 3.0  # only the continuous regressor has powers
    z = np.column_stack([np.ones(y.size), x[:, 0], dummy])
    expected, _, _ = nested_f(y, z, np.column_stack([x[:, 0] ** 2, x[:, 0] ** 3, x[:, 0] ** 4]))
    assert result.statistic == pytest.approx(expected, rel=1e-8)


def test_the_test_rejects_a_missing_square_and_not_a_correct_model() -> None:
    y, x = data(seed=1, curvature=0.8)
    assert reset(y, x).rejects(0.001)
    y, x = data(seed=2)
    assert not reset(y, x).rejects(0.001)


def test_a_collinear_regressor_changes_nothing_and_collinear_powers_are_left_out() -> None:
    y, x = data(curvature=0.4)
    same = reset(y, np.column_stack([x, 2 * x[:, 0] + x[:, 1]]))
    assert same.statistic == pytest.approx(reset(y, x).statistic, rel=1e-8)
    duplicated = reset(y, np.column_stack([x[:, 0], x[:, 0]]), rhs=True)
    assert duplicated.df[0] == 3.0  # the 6 powers of two identical columns span only 3
    assert "3 of the 6 powers" in duplicated.summary("en")
    dummy = (x[:, 0] > 0).astype(float)
    for rhs in (False, True):  # the powers of an indicator are the indicator itself
        with pytest.raises(NumericalError, match="collinear"):
            reset(y, dummy, rhs=rhs)


def test_the_statistic_does_not_depend_on_the_scale_of_the_regressors() -> None:
    y, x = data(curvature=0.4)
    base = reset(y, x)
    assert reset(y, x * np.array([1000.0, 0.001])).statistic == pytest.approx(
        base.statistic, rel=1e-6
    )
    assert reset(y, x, rhs=True).statistic == pytest.approx(
        reset(y, x * np.array([1000.0, 0.001]), rhs=True).statistic, rel=1e-6
    )


def test_shapes_and_a_perfect_fit_are_refused() -> None:
    y, x = data()
    with pytest.raises(SpecificationError, match="one row per observation"):
        reset(y, x[:-1])
    with pytest.raises(NumericalError, match="perfectly"):
        reset(x[:, 0] ** 2 + x[:, 0], x[:, 0])


def test_by_name_and_alias_in_both_languages() -> None:
    y, x = data()
    result = run_test("ovtest", y, x)
    assert result.name == "reset"
    assert "Ramsey RESET test for omitted variables" in result.summary("en")
    assert "Test RESET de Ramsey" in result.summary("fr")
