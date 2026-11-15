"""Correlation matrix between units, and the tolerance the tests do not share."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel.diagnostics import (
    cross_correlation,
    lm_cross_section,
    run_test,
)


def declare(wide: np.ndarray) -> xt.PanelData:
    """Declare a ``(periods, units)`` array (NaN = absent) as a panel."""
    periods, units = wide.shape
    frame = pd.DataFrame(
        {
            "id": np.tile(np.arange(1, units + 1), periods),
            "t": np.repeat(np.arange(1, periods + 1), units),
            "v": wide.ravel(),
        }
    )
    return xt.xtset(frame.dropna(subset=["v"]).reset_index(drop=True), "id", "t")


def factor_data(seed: int, periods: int = 12, units: int = 8, load: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.normal(size=(periods, units)) + load * rng.normal(size=(periods, 1))


def test_the_matrix_equals_numpy_on_a_balanced_panel() -> None:
    wide = factor_data(8)
    result = cross_correlation(declare(wide), "v")
    assert np.allclose(result.correlation.to_numpy(), np.corrcoef(wide, rowvar=False), atol=1e-12)
    assert (np.diag(result.overlap.to_numpy()) == 12).all()
    rho = result.correlation.to_numpy()[np.triu_indices(8, 1)]
    assert result.average == pytest.approx(rho.mean()) and result.n_pairs == 28
    assert result.average_abs == pytest.approx(np.abs(rho).mean()) and result.n_undefined == 0
    assert result.to_frame() is result.correlation


def test_short_or_constant_pairs_are_left_out_not_refused() -> None:
    wide = factor_data(9, periods=8, units=5)
    wide[:6, 0] = np.nan  # unit 1 has 2 periods: its 4 pairs have 2 common periods
    wide[:, 4] = 3.0  # unit 5 is constant: 3 more pairs have no variation
    result = cross_correlation(declare(wide), "v")
    assert result.n_undefined == 4 + 3 and result.n_pairs == 3
    assert result.correlation.iloc[0, 1:].isna().all()
    assert "7 pair(s)" in result.summary("en") and "7 paire(s)" in result.summary("fr")
    with pytest.raises(SpecificationError, match="min_overlap must be at least 2"):
        cross_correlation(declare(wide), "v", min_overlap=1)


def test_the_tests_refuse_what_the_matrix_tolerates() -> None:
    wide = factor_data(10, periods=8, units=5)
    wide[:, 4] = 3.0
    with pytest.raises(NumericalError, match="no variation over their common periods"):
        lm_cross_section(declare(wide), "v")


def test_by_name_and_alias_in_both_languages() -> None:
    panel = declare(factor_data(11))
    assert run_test("xttest2", panel, "v").statistic == lm_cross_section(panel, "v").statistic
    assert "Breusch-Pagan LM test" in run_test("bp_lm", panel, "v").summary("en")
    assert "Test de Friedman" in run_test("friedman", panel, "v").summary("fr")
    assert "Test de Frees" in run_test("frees", panel, "v").summary("fr")
