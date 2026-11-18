"""A figure called wrongly is refused with a message that says what the figure accepts (M.2).

Without the check, ``panel.plot("states")`` ended in Python's bare ``TypeError``, which names
neither the figure nor the arguments it takes.  The same check as the one of ``run_panel_test``
runs before any drawing, for every registered kind.
"""

from __future__ import annotations

import inspect
import re

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")

import pyxtpanel as xt  # noqa: E402
from pyxtpanel._errors import SpecificationError  # noqa: E402
from pyxtpanel.plot import KINDS  # noqa: E402


@pytest.fixture(scope="module")
def panel() -> xt.PanelData:
    rng = np.random.default_rng(1)
    frame = pd.DataFrame(
        {"id": np.repeat(np.arange(1, 6), 6), "t": np.tile(np.arange(2000, 2006), 5)}
    )
    frame["x"] = rng.normal(size=len(frame))
    return xt.xtset(frame, "id", "t")


def _required(kind: str) -> list[str]:
    parameters = list(inspect.signature(KINDS.get(kind)).parameters.values())[1:]
    return [p.name for p in parameters if p.default is inspect.Parameter.empty]


@pytest.mark.parametrize("kind", KINDS.names())
def test_an_unknown_option_is_refused_with_the_arguments_of_the_figure(
    panel: xt.PanelData, kind: str
) -> None:
    given = {name: "x" for name in _required(kind)}
    with pytest.raises(SpecificationError) as error:
        panel.plot(kind, bogus_option=1, **given)
    message = str(error.value)
    assert f"figure {kind!r}" in message and "bogus_option" in message
    assert "Arguments of the figure:" in message and "lang=" in message
    assert getattr(error.value, "key", None) == "plot.bad_call"


@pytest.mark.parametrize("kind", [k for k in KINDS.names() if _required(k)])
def test_a_missing_required_argument_is_named(panel: xt.PanelData, kind: str) -> None:
    with pytest.raises(SpecificationError) as error:
        panel.plot(kind)
    message = str(error.value)
    assert f"figure {kind!r}" in message
    # Python words it "keyword-only" for an option that can only be given by name
    name = _required(kind)[0]
    assert re.search(rf"missing a required (keyword-only )?argument: '{name}'", message)


def test_the_states_figure_without_a_variable_is_the_reported_case(panel: xt.PanelData) -> None:
    with pytest.raises(SpecificationError, match="variable"):
        panel.plot("states")


def test_the_message_follows_the_language_asked_for_the_figure(panel: xt.PanelData) -> None:
    with pytest.raises(SpecificationError, match="Arguments de la figure") as french:
        panel.plot("participation", lang="fr", bogus_option=1)
    assert "Arguments of the figure" not in str(french.value)
    with pytest.raises(SpecificationError, match="Arguments of the figure"):
        panel.plot("participation", lang="en", bogus_option=1)


def test_a_valid_call_still_draws(panel: xt.PanelData) -> None:
    figure = panel.plot("participation")
    assert figure.axes
