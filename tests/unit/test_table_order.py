"""The rows of the small result tables keep the order the package writes them in, on every pandas.

pandas before 2.2 sorts the index of ``DataFrame({"col": a_dict})`` alphabetically; the figures of
``unbalancedness``, the key figures of the audit and the quantiles of ``complete_cases`` are built
from ordered dictionaries, so a reordering would change what a user reads (found by running the
suite on pandas 2.1 and numpy 1.26, J.1.3).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import pyxtpanel as xt


def _panel() -> xt.PanelData:
    rng = np.random.default_rng(1)
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(6), 5),
            "t": np.tile(np.arange(5), 6),
            "x": rng.normal(size=30),
            "y": rng.normal(size=30),
        }
    ).drop(index=[3, 12, 13])
    return xt.xtset(frame, "id", "t")


def test_unbalancedness_figures_are_in_the_documented_order() -> None:
    index = list(_panel().unbalancedness().to_frame().index)
    assert index == ["gamma", "nu", "T_mean", "T_harmonic", "T_min", "T_max", "fill_rate"]


def test_audit_key_figures_start_with_the_sizes_and_end_with_the_row_share() -> None:
    facts = list(_panel().audit().facts.index)
    assert facts[:4] == ["n_units", "n_periods", "n_obs", "T_min"]
    assert facts[-2:] == ["units_single_obs", "complete_row_share"]


def test_complete_cases_quantiles_run_from_the_mean_to_the_maximum() -> None:
    index = list(_panel().complete_cases(["x", "y"]).ti.index)
    assert index[:2] == ["mean", "min"] and index[-1] == "max"
