"""Sargan test: values computed by hand and the role of the scale."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics import run_test
from pyxtpanel.diagnostics.postestimation.hansen import hansen
from pyxtpanel.diagnostics.postestimation.sargan import sargan


def test_hand_computed_with_a_scale() -> None:
    # m'Wm = 5 with W = I, divided by the scale 2: 2.5; chi2(1) tail = erfc(sqrt(2.5 / 2)).
    result = sargan([1.0, 2.0], np.eye(2), 1, scale=2.0)
    assert result.statistic == 2.5 and result.df == (1.0,) and result.name == "sargan"
    assert result.pvalue == pytest.approx(math.erfc(math.sqrt(1.25)), rel=1e-12)
    assert result.options["scale"] == 2.0


def test_with_unit_scale_it_is_the_hansen_quadratic_form() -> None:
    weight = np.array([[2.0, 1.0], [1.0, 3.0]])
    assert (
        sargan([1.0, -1.0], weight, 1, scale=1.0).statistic
        == hansen([1.0, -1.0], weight, 1).statistic
    )


def test_the_scale_is_required_and_must_be_positive() -> None:
    with pytest.raises(TypeError):
        sargan([1.0, 2.0], np.eye(2), 1)  # type: ignore[call-arg]
    for bad in (0.0, -1.0, math.nan, math.inf):
        with pytest.raises(SpecificationError, match="scale must be positive"):
            sargan([1.0, 2.0], np.eye(2), 1, scale=bad)


def test_not_over_identified_gives_nan() -> None:
    result = sargan([1.0], np.eye(1), 1, scale=1.0)
    assert math.isnan(result.statistic) and math.isnan(result.pvalue)


def test_two_stage_least_squares_by_hand() -> None:
    # Sargan after 2SLS is N * R^2 of the residuals on the instruments: check on a tiny design.
    rng = np.random.default_rng(11)
    n = 60
    z = np.column_stack([rng.normal(size=(n, 3)), np.ones(n)])
    x = z[:, :3] @ np.array([0.5, 0.4, 0.3]) + rng.normal(size=n)
    y = 1.0 + 0.8 * x + rng.normal(size=n)
    xm = np.column_stack([x, np.ones(n)])
    w = np.linalg.inv(z.T @ z)
    beta = np.linalg.solve(xm.T @ z @ w @ z.T @ xm, xm.T @ z @ w @ z.T @ y)
    e = y - xm @ beta
    fitted = z @ np.linalg.solve(z.T @ z, z.T @ e)
    n_r2 = n * float(fitted @ fitted) / float(e @ e)
    result = run_test("sargan", z.T @ e, w, 2, scale=float(e @ e) / n)
    assert result.statistic == pytest.approx(n_r2, rel=1e-10)
    assert result.df == (2.0,)
