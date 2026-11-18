"""The level of agreement that closed forms can reach, measured against Stata itself (B.3.4).

``validation/self_agreement.py`` ran the battery of the campaign in Stata on a twin of the frozen
data (units renamed in another order, rows shuffled) and froze three measures in
``stata_self_agreement.json``.  These tests read the file (no Stata) and recompute the one that
needs only the package, so that the tolerance the validation claims cannot drift from the data.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pytest
from validation import campaign
from validation.self_agreement import error

GOLDEN = Path(__file__).resolve().parent
MEASURE = json.loads((GOLDEN / "stata_self_agreement.json").read_text(encoding="utf-8"))
BATTERY = np.array(
    json.loads((GOLDEN / "campaign_golden.json").read_text(encoding="utf-8"))["cases"]["campaign"][
        "battery"
    ],
    dtype=float,
)


def test_stata_agrees_with_itself_far_below_the_tolerance_of_the_closed_forms() -> None:
    """Renaming the units changes the order of Stata's additions: the noise is ~1e-14."""
    assert MEASURE["stata_vs_twin"]["max"] < 1e-12
    assert MEASURE["stata_vs_twin"]["cells"] > 10_000


def test_the_package_agrees_with_stata_to_about_1e_12() -> None:
    assert MEASURE["package_vs_stata"]["max"] < 1e-11
    assert MEASURE["package_vs_stata"]["cells_above_1e-12"] == 0


def test_reading_a_value_on_the_screen_loses_far_more_than_the_package_does() -> None:
    """The log shows rounded values: the frozen references are taken from ``r()`` instead."""
    package = MEASURE["package_vs_stata"]["max"]
    assert MEASURE["display"]["%12.0g"]["max"] > 1000 * package
    assert MEASURE["display"]["%9.0g"]["max"] > MEASURE["display"]["%12.0g"]["max"]


@pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")
def test_the_measure_of_the_package_is_still_true_today() -> None:
    """Recompute the package side (Stata's side is frozen): the claim must keep holding."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        package = np.vstack([campaign.python_row(s.frame, s.delta, s.fmt) for s in campaign.load()])
    errors = error(package, BATTERY)
    assert np.isfinite(errors).all()
    assert errors.max() < 1e-11
    assert errors.max() == pytest.approx(MEASURE["package_vs_stata"]["max"], abs=1e-11)
