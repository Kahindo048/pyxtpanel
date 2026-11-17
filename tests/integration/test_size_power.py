"""Size and power of the residual and specification tests by simulation (plan L.8.12).

The scenarios are those of ``validation/simulate.py``, run here with few replications: the bounds
are wide enough for the Monte Carlo error (about 1.3 points at 300 replications) and narrow
enough to catch a statistic or a p-value that is wrong.  A test known to be oversized
(``wald_groupwise``, as Stata's ``xttest3``) has bounds that say so.
"""

from __future__ import annotations

import pytest
from validation import simulate

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


def test_the_panel_durbin_watson_statistic_is_about_two_without_autocorrelation() -> None:
    null, alternative, (low, high), upper = simulate.DURBIN_WATSON
    assert low <= simulate.mean_statistic(null, 100, 7) <= high
    assert simulate.mean_statistic(alternative, 100, 8) <= upper
