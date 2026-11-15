"""Variance inflation factors: the diagonal of the inverse correlation matrix as oracle."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import VarianceInflation, vif


def regressors(n: int = 500) -> np.ndarray:
    rng = np.random.default_rng(0)
    a, b = rng.normal(size=(2, n))
    return np.column_stack([a, 0.8 * a + 0.6 * b, rng.normal(size=n) + 0.5 * a])


def test_the_factors_are_the_diagonal_of_the_inverse_correlation_matrix() -> None:
    x = regressors()
    expected = np.diag(np.linalg.inv(np.corrcoef(x, rowvar=False)))
    result = vif(x, sort=False)
    assert result.vif.to_numpy() == pytest.approx(expected, rel=1e-10)
    assert result.mean == pytest.approx(expected.mean(), rel=1e-12)


def test_the_uncentered_factors_include_the_constant() -> None:
    x = np.column_stack([regressors() + 3.0])
    z = np.column_stack([x, np.ones(x.shape[0])])
    scaled = z / np.sqrt(np.sum(z**2, axis=0))
    expected = np.diag(np.linalg.inv(scaled.T @ scaled))
    result = vif(x, uncentered=True, sort=False)
    assert result.vif.to_numpy() == pytest.approx(expected, rel=1e-10)
    assert list(result.vif.index) == ["x1", "x2", "x3", "_cons"]


def test_sorting_names_and_the_table_form() -> None:
    frame = pd.DataFrame(regressors(), columns=["a", "b", "c"])
    result = vif(frame)
    assert isinstance(result, VarianceInflation)
    assert list(result.vif) == sorted(result.vif, reverse=True)
    assert set(result.vif.index) == {"a", "b", "c"}
    table = result.to_frame()
    assert list(table.columns) == ["VIF", "1/VIF"]
    assert table["1/VIF"].to_numpy() == pytest.approx(1 / table["VIF"].to_numpy())
    assert list(vif(frame, sort=False).vif.index) == ["a", "b", "c"]
    assert list(vif(regressors(), names=["p", "q", "r"], sort=False).vif.index) == ["p", "q", "r"]


def test_independent_regressors_have_a_factor_close_to_one() -> None:
    x = np.random.default_rng(1).normal(size=(5000, 3))
    assert vif(x).vif.max() < 1.01


def test_summary_in_both_languages() -> None:
    result = vif(pd.DataFrame(regressors(), columns=["a", "b", "c"]))
    assert "Variance inflation factors (centered): 3 variables" in result.summary("en")
    assert "Facteurs d'inflation de la variance (centrés)" in result.summary("fr")
    assert "(uncentered)" in vif(regressors(), uncentered=True).summary("en")


def test_refusals() -> None:
    x = regressors()
    with pytest.raises(SpecificationError, match="at least 2"):
        vif(x[:, :1])
    with pytest.raises(SpecificationError, match="one entry per column"):
        vif(x, names=["a"])
    with pytest.raises(NumericalError, match="constant"):
        vif(np.column_stack([x, np.ones(x.shape[0])]))
    with pytest.raises(NumericalError, match="linear combination"):
        vif(np.column_stack([x, x[:, 0] + x[:, 1]]))
