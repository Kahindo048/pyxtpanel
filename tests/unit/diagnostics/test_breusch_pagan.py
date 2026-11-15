"""Breusch-Pagan / Cook-Weisberg test: hand-computed values, identities and refusals."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import breusch_pagan, run_test
from pyxtpanel.diagnostics._shared.ols import ols_fit


def data(seed: int = 5, n: int = 80) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    z = rng.normal(size=n)
    return rng.normal(size=n) * np.exp(0.4 * z), z


def test_hand_computed_perfect_dependence() -> None:
    # e^2 = (1, 1, 9, 9), mean 5: u = (0.2, 0.2, 1.8, 1.8) = 0.2 + 1.6 z for z = (0, 0, 1, 1).
    # TSS = 4 * 0.8^2 = 2.56 = MSS, R^2 = 1: normal = 1.28, iid = N R^2 = 4, both on 1 df.
    e, z = np.array([1.0, -1.0, 3.0, -3.0]), np.array([0.0, 0.0, 1.0, 1.0])
    normal = breusch_pagan(e, z)
    assert normal.statistic == pytest.approx(1.28, rel=1e-13) and normal.df == (1.0,)
    assert normal.pvalue == pytest.approx(math.erfc(math.sqrt(0.64)), rel=1e-12)
    iid = breusch_pagan(e, z, variant="iid")
    assert iid.statistic == pytest.approx(4.0, rel=1e-13)
    assert iid.pvalue == pytest.approx(math.erfc(math.sqrt(2.0)), rel=1e-12)
    assert normal.statistic_name == "chi2" and normal.n_obs == 4


def test_the_three_forms_follow_from_the_r_squared_of_a_simple_regression() -> None:
    e, z = data()
    u = e**2 / np.mean(e**2)
    r2 = float(np.corrcoef(u, z)[0, 1] ** 2)
    n = e.size
    tss = float(np.sum((u - u.mean()) ** 2))
    assert breusch_pagan(e, z).statistic == pytest.approx(r2 * tss / 2, rel=1e-11)
    assert breusch_pagan(e, z, variant="iid").statistic == pytest.approx(n * r2, rel=1e-11)
    f = breusch_pagan(e, z, variant="fstat")
    assert f.statistic == pytest.approx(r2 / (1 - r2) * (n - 2), rel=1e-10)
    assert f.df == (1.0, float(n - 2)) and f.statistic_name == "F"


def test_several_variables_count_their_rank_not_their_number() -> None:
    e, z = data()
    both = np.column_stack([z, z, np.ones_like(z), z**2])  # a copy, a constant and a square
    result = breusch_pagan(e, both)
    assert result.df == (2.0,)
    assert result.statistic == pytest.approx(
        breusch_pagan(e, np.column_stack([z, z**2])).statistic, rel=1e-10
    )


def test_invariant_to_the_scale_of_the_residuals_and_the_row_order() -> None:
    e, z = data()
    base = breusch_pagan(e, z)
    assert breusch_pagan(7.5 * e, z).statistic == pytest.approx(base.statistic, rel=1e-12)
    order = np.random.default_rng(1).permutation(e.size)
    assert breusch_pagan(e[order], z[order]).statistic == pytest.approx(base.statistic, rel=1e-12)


def test_a_variance_that_depends_on_the_variable_is_rejected() -> None:
    e, z = data(n=500)
    assert breusch_pagan(e, z).rejects(0.01)
    assert not breusch_pagan(np.random.default_rng(2).normal(size=500), z).rejects(0.01)


def test_refusals() -> None:
    e, z = data()
    with pytest.raises(SpecificationError, match="variant must be 'normal', 'iid' or 'fstat'"):
        breusch_pagan(e, z, variant="robust")  # type: ignore[arg-type]
    with pytest.raises(SpecificationError, match="one row per residual"):
        breusch_pagan(e, z[:-1])
    with pytest.raises(NumericalError, match="all residuals are zero"):
        breusch_pagan(np.zeros(5), np.arange(5.0))
    with pytest.raises(NumericalError, match="constant or collinear"):
        breusch_pagan(e, np.ones_like(e))
    with pytest.raises(NumericalError, match="too few"):
        breusch_pagan(
            np.array([1.0, 2.0, 3.0]),
            np.array([[0.0, 1.0], [1.0, 0.0], [2.0, 5.0]]),
            variant="fstat",
        )


def test_by_name_in_both_languages() -> None:
    e, z = data()
    result = run_test("breusch_pagan", e, z)
    assert "Breusch-Pagan / Cook-Weisberg" in result.summary("en")
    assert "variance constante" in result.summary("fr")


def test_the_auxiliary_regression_drops_constant_and_collinear_columns() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(30, 2))
    y = x @ np.array([1.0, 2.0]) + rng.normal(size=30)
    fit = ols_fit(y, np.column_stack([x, x[:, 0] + x[:, 1], np.full(30, 4.0)]))
    assert fit.df_model == 2 and fit.df_resid == 27
    coefficients, *_ = np.linalg.lstsq(np.column_stack([np.ones(30), x]), y, rcond=None)
    rss = float(np.sum((y - np.column_stack([np.ones(30), x]) @ coefficients) ** 2))
    assert fit.rss == pytest.approx(rss, rel=1e-12)
    assert fit.mss + fit.rss == pytest.approx(fit.tss, rel=1e-12)
    assert fit.r2 == pytest.approx(fit.mss / fit.tss) and fit.f_statistic > 0
    flat = ols_fit(y, np.ones((30, 1)))
    assert flat.df_model == 0 and flat.mss == 0.0 and math.isnan(flat.f_statistic)
