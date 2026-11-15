"""White's test: auxiliary-regression oracle, collinear terms, invariances and refusals."""

from __future__ import annotations

import numpy as np
import pytest

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import run_test, white


def oracle(e: np.ndarray, x: np.ndarray) -> tuple[float, int]:
    """N R^2 and the rank of the auxiliary regression (least squares, explicit constant)."""
    x = x.reshape(e.size, -1)
    centred = x - x.mean(axis=0)
    terms = [centred[:, i] * centred[:, j] for i in range(x.shape[1]) for j in range(i + 1)]
    design = np.column_stack([np.ones(e.size), x, *terms])
    y = e**2
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    rss = float(np.sum((y - design @ beta) ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    return e.size * (1 - rss / tss), int(np.linalg.matrix_rank(design)) - 1


def data(p: int, seed: int = 4, n: int = 200) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, p))
    return rng.normal(size=n) * np.exp(0.5 * x[:, 0]), x


@pytest.mark.parametrize("p", [1, 2, 3, 4])
def test_it_is_n_r_squared_of_the_auxiliary_regression(p: int) -> None:
    e, x = data(p)
    statistic, rank = oracle(e, x)
    result = white(e, x)
    assert result.statistic == pytest.approx(statistic, rel=1e-10)
    assert result.df == (float(rank),) and rank == p + p * (p + 1) // 2
    assert result.n_obs == 200 and result.statistic_name == "chi2"


def test_a_dummy_regressor_has_a_collinear_square_that_is_dropped() -> None:
    rng = np.random.default_rng(8)
    x = np.column_stack([rng.normal(size=300), (rng.random(300) < 0.4).astype(float)])
    e = rng.normal(size=300) * np.exp(0.5 * x[:, 0])
    # 2 regressors: x1, d, x1^2, x1*d, d^2 -> d^2 is a copy of d once centred terms are used.
    assert white(e, x).df == (4.0,)
    assert white(e, x).statistic == pytest.approx(oracle(e, x)[0], rel=1e-9)


def test_invariant_to_the_scale_of_the_residuals_and_to_the_row_order() -> None:
    e, x = data(2)
    base = white(e, x)
    assert white(3.0 * e, x).statistic == pytest.approx(base.statistic, rel=1e-11)
    order = np.random.default_rng(3).permutation(e.size)
    assert white(e[order], x[order]).statistic == pytest.approx(base.statistic, rel=1e-11)


def test_rejects_heteroskedasticity_and_not_homoskedasticity() -> None:
    e, x = data(1, n=600)
    assert white(e, x).rejects(0.01)
    assert not white(np.random.default_rng(7).normal(size=600), x).rejects(0.01)


def test_refusals() -> None:
    e, x = data(1)
    with pytest.raises(SpecificationError, match="one row per residual"):
        white(e, x[:-1])
    with pytest.raises(SpecificationError, match="at least one regressor"):
        white(e, np.empty((e.size, 0)))
    with pytest.raises(NumericalError, match="constant or collinear"):
        white(e, np.ones((e.size, 1)))


def test_by_name_in_both_languages() -> None:
    e, x = data(1)
    result = run_test("white", e, x)
    assert "White's test for heteroskedasticity" in result.summary("en")
    assert "homoscédasticité" in result.summary("fr")
