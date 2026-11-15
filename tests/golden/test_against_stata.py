"""Every frozen Stata result against the package (no Stata needed).

The test is generic: it reads each ``tests/golden/<domain>_golden.json`` written by
``python -m validation.run <domain>`` and runs the Python side of every registered case of that
domain.  A new validation case needs no new test.
"""

import json
from pathlib import Path

import numpy as np
import pytest
from validation import compare
from validation.registry import CASES, discover

pytestmark = pytest.mark.golden

GOLDEN_DIR = Path(__file__).parent
MASTERS = {
    path.name.removesuffix("_golden.json"): json.loads(path.read_text(encoding="utf-8"))
    for path in sorted(GOLDEN_DIR.glob("*_golden.json"))
}
discover()
PARAMS = [(domain, tag) for domain, master in MASTERS.items() for tag in master["cases"]]


@pytest.mark.parametrize(("domain", "tag"), PARAMS, ids=[f"{d}:{t}" for d, t in PARAMS])
def test_case_agrees_with_stata(domain: str, tag: str) -> None:
    case = CASES[tag]
    assert case.domain == domain
    stata = {k: np.array(v, dtype=float) for k, v in MASTERS[domain]["cases"][tag].items()}
    python = {k: np.asarray(v, dtype=float) for k, v in case.python().items()}
    diffs = compare.compare_case(python, stata, case.atol, case.rtol)
    assert all(d.ok for d in diffs), compare.report(tag, diffs)


@pytest.mark.parametrize("domain", sorted(MASTERS))
def test_every_registered_case_is_frozen(domain: str) -> None:
    registered = {tag for tag, case in CASES.items() if case.domain == domain}
    assert registered == set(MASTERS[domain]["cases"]), "re-run python -m validation.run"


def test_the_golden_files_say_which_stata_produced_them() -> None:
    for master in MASTERS.values():
        assert master["meta"]["stata"].startswith("Stata 17")
