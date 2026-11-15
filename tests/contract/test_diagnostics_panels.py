"""What each test on a panel declares about balance is what it does (L.1.5), and aliases are clean.

A test declared ``balanced="required"`` must refuse an unbalanced panel with the uniform message
that names the first incomplete unit; a test declared ``"both"`` must accept it.  Every category A
test needs an entry in ``RECIPES`` (the call that runs it on the common panel), so a new test
cannot be added without being exercised here.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

import pyxtpanel as xt
from pyxtpanel._errors import DataError
from pyxtpanel.diagnostics import TestResult, list_tests, run_test
from pyxtpanel.diagnostics._core import TESTS

#: test name -> (positional arguments after the panel, options); variables are ``y`` and ``x``.
RECIPES: dict[str, tuple[tuple[Any, ...], dict[str, Any]]] = {
    "pesaran_cd": (("y",), {}),
    "wooldridge": (("y", ["x"]), {}),
    "harris_tzavalis": (("y",), {}),
    "levin_lin_chu": (("y",), {"lags": 1}),
    "im_pesaran_shin": (("y",), {}),
    "fisher": (("y",), {"lags": 1}),
    "lm_cross_section": (("y",), {}),
    "friedman": (("y",), {}),
    "frees": (("y",), {}),
    "slope_homogeneity": (("y", ["x"]), {}),
}


def make_panel(balanced: bool) -> xt.PanelData:
    rng = np.random.default_rng(0)
    n, t = 30, 24
    y = rng.normal(size=(n, t)).cumsum(axis=1)
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, n + 1), t),
            "t": np.tile(np.arange(1, t + 1), n),
            "y": y.ravel(),
            "x": rng.normal(size=n * t),
        }
    )
    if not balanced:  # units 1 to 10 leave early or enter late, never with a hole
        late = (frame["id"] <= 5) & (frame["t"] <= 3)
        early = (frame["id"].between(6, 10)) & (frame["t"] > 20)
        frame = frame[~(late | early)].reset_index(drop=True)
    return xt.xtset(frame, "id", "t")


def run(name: str, panel: xt.PanelData) -> TestResult:
    args, options = RECIPES[name]
    return run_test(name, panel, *args, **options)


CATEGORY_A = [name for name, spec in TESTS.specs() if spec.category == "A"]


def test_every_test_on_a_panel_has_a_recipe() -> None:
    assert sorted(CATEGORY_A) == sorted(RECIPES)


@pytest.mark.parametrize("name", CATEGORY_A)
def test_a_balanced_panel_is_accepted_by_every_test(name: str) -> None:
    assert np.isfinite(run(name, make_panel(True)).statistic)


@pytest.mark.parametrize("name", CATEGORY_A)
def test_the_balance_a_test_declares_is_the_balance_it_enforces(name: str) -> None:
    spec = TESTS.spec(name)
    panel = make_panel(False)
    if spec.balanced == "required":
        with pytest.raises(DataError, match="needs a strongly balanced panel: unit 1 is observed"):
            run(name, panel)
    else:
        assert spec.balanced == "both"
        assert np.isfinite(run(name, panel).statistic)


def test_the_listing_says_which_tests_need_a_balanced_panel() -> None:
    listing = list_tests()
    required = set(listing.index[listing["balanced"] == "required"])
    assert required == {"harris_tzavalis", "levin_lin_chu"}


def test_aliases_are_lower_case_and_never_collide() -> None:
    names = set(TESTS.names())
    seen: dict[str, str] = {}
    for name, spec in TESTS.specs():
        for alias in spec.aliases:
            assert alias == alias.lower() and alias not in names, (name, alias)
            assert alias not in seen, f"{alias!r} is an alias of {seen[alias]} and of {name}"
            seen[alias] = name


def test_stata_subcommand_names_are_aliases() -> None:
    listing = list_tests()
    for command, name in (
        ("xtunitroot ht", "harris_tzavalis"),
        ("xtunitroot llc", "levin_lin_chu"),
        ("xtunitroot ips", "im_pesaran_shin"),
    ):
        assert listing.loc[name, "stata"] == command
        assert command.split()[-1] in TESTS.spec(name).aliases
