"""``usable_equations`` agrees with the observations ``pyxtabond2`` uses (K.2.4, K.7.1).

The numbers come from ``python -m validation.cross_pyxtabond2``: one lag of the dependent
variable, first differences (instruments from lag 2) and orthogonal deviations (from lag 1).
"""

import json
from pathlib import Path

import pytest
from validation import panels

import pyxtpanel as xt

FROZEN = json.loads((Path(__file__).parent / "equations_pyxtabond2.json").read_text("utf-8"))


@pytest.mark.parametrize("transform", ["fd", "fod"])
@pytest.mark.parametrize("name", sorted(FROZEN))
def test_the_equations_counted_are_the_observations_pyxtabond2_uses(
    name: str, transform: str
) -> None:
    panel = xt.xtset(panels.load(name)[["id", "t", "x"]], "id", "t")
    result = panel.usable_equations("x", lags=1, transform=transform)  # type: ignore[arg-type]
    assert result.n_instrumented == FROZEN[name][transform]
    assert result.n_usable == FROZEN[name][transform]


def test_first_differences_and_orthogonal_deviations_differ_where_there_are_gaps() -> None:
    assert FROZEN["gaps"]["fod"] > FROZEN["gaps"]["fd"]
    assert FROZEN["balanced"]["fod"] == FROZEN["balanced"]["fd"]
