"""Size and power of the panel and residual tests by simulation (plans L.7.5 and L.8.12).

The scenarios are those of ``validation/simulate.py``, run here with few replications: the bounds
are wide enough for the Monte Carlo error (about 1.3 points at 300 replications) and narrow
enough to catch a statistic or a p-value that is wrong.  A test known to be oversized
(``wald_groupwise``, as Stata's ``xttest3``; ``levin_lin_chu`` and ``harris_tzavalis``, as
``xtunitroot``) has bounds that say so.
"""

from __future__ import annotations

import pytest
from validation import simulate

from pyxtpanel.diagnostics._core.spec import list_tests

REPS = 300

pytestmark = pytest.mark.slow


@pytest.mark.parametrize("name", list(simulate.STUDIES))
def test_the_rejection_rate_under_the_null_is_near_the_level(name: str) -> None:
    study = simulate.STUDIES[name]
    seed = 1000 + list(simulate.STUDIES).index(name)
    size = simulate.rejection_rate(study.null, REPS, seed)
    low, high = study.size_bounds
    assert low <= size <= high, f"{name}: size {size:.3f} outside [{low}, {high}] ({study.note})"


@pytest.mark.parametrize("name", list(simulate.STUDIES))
def test_the_rejection_rate_under_the_alternative_is_high(name: str) -> None:
    study = simulate.STUDIES[name]
    seed = 1500 + list(simulate.STUDIES).index(name)
    power = simulate.rejection_rate(study.alternative, REPS, seed)
    assert power >= study.min_power, f"{name}: power {power:.3f} below {study.min_power}"


def test_every_test_of_category_a_has_a_size_and_power_study() -> None:
    # Plan L.7.5: a category A test added without its simulation fails here.  The Durbin-Watson
    # statistic has no p-value and its own study below.
    table = list_tests()
    category_a = set(table.index[table["category"] == "A"])
    missing = category_a - set(simulate.STUDIES) - {"panel_durbin_watson"}
    assert not missing, f"no size and power study for {sorted(missing)}"


def test_the_panel_durbin_watson_statistic_is_about_two_without_autocorrelation() -> None:
    null, alternative, (low, high), upper = simulate.DURBIN_WATSON
    assert low <= simulate.mean_statistic(null, 100, 7) <= high
    assert simulate.mean_statistic(alternative, 100, 8) <= upper
