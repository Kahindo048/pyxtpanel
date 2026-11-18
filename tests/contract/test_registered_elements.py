"""Every registered element is complete: tested, summarised, and validated as it claims (A.5.10).

The registries (operators, descriptions, exploration tools, audit checks, figures) let a feature
be added in one file.  These meta-tests fail when a new element arrives without what the others
have: a test that names it, a one-line summary, and an explicit statement of how it is validated.
The statistical tests have their own gate (``test_diagnostics_tests.py``), the descriptions and
operators their provenance checks (``test_public_api.py``), the tools their validation claims
(``test_explore_tools.py``).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from pyxtpanel import _i18n
from pyxtpanel.describe._spec import DESCRIPTIONS
from pyxtpanel.explore._spec import TOOLS
from pyxtpanel.explore.checks._spec import CHECKS
from pyxtpanel.ops._spec import OPERATORS
from pyxtpanel.plot import KINDS

ROOT = Path(__file__).resolve().parents[2]
TESTS = {
    path: path.read_text(encoding="utf-8", errors="ignore")
    for path in (ROOT / "tests").rglob("*.py")
    if path.name != Path(__file__).name
}

REGISTRIES = {
    "operator": OPERATORS,
    "description": DESCRIPTIONS,
    "tool": TOOLS,
    "check": CHECKS,
    "figure": KINDS,
}
ELEMENTS = [(label, name) for label, registry in REGISTRIES.items() for name in registry.names()]

#: Operators whose results were compared with Stata 17 (``operators`` and ``transforms`` golden
#: domains and the random campaign).  Any other operator is validated by values worked out by hand
#: and against pandas (``tests/properties/test_cross_pandas.py``); its ``stata`` field then names
#: the command it is the equivalent of, not a comparison that was run.
OPERATORS_COMPARED_WITH_STATA = {
    "lag",
    "lead",
    "difference",
    "seasonal",
    "within",
    "between",
    "fod",
    "quasi_demean",
    "twoway",
}


def _named_in_a_test(name: str) -> list[Path]:
    pattern = re.compile(rf"""["'(.\s]{re.escape(name)}["'(\s]""")
    return [path for path, text in TESTS.items() if pattern.search(text)]


@pytest.mark.parametrize(("label", "name"), ELEMENTS, ids=lambda v: str(v))
def test_every_registered_element_is_named_in_a_test(label: str, name: str) -> None:
    assert _named_in_a_test(name), f"no test mentions the {label} {name!r}"


@pytest.mark.parametrize(("label", "name"), ELEMENTS, ids=lambda v: str(v))
def test_every_registered_element_has_a_summary(label: str, name: str) -> None:
    summary = REGISTRIES[label].spec(name).summary
    assert summary.strip() and not summary.startswith("TODO"), f"{label} {name!r}: no summary"


def test_every_operator_says_how_it_is_validated() -> None:
    """A new operator must be classified: compared with Stata, or validated by hand and pandas."""
    names = set(OPERATORS.names())
    assert names >= OPERATORS_COMPARED_WITH_STATA, "a Stata-validated operator was removed"
    for name in sorted(names - OPERATORS_COMPARED_WITH_STATA):
        assert _named_in_a_test(name), (
            f"operator {name!r} is neither compared with Stata (add it to "
            "OPERATORS_COMPARED_WITH_STATA once a case exists) nor tested by hand"
        )


@pytest.mark.parametrize("name", sorted(set(DESCRIPTIONS.names())))
def test_a_description_that_claims_stata_has_a_validation_case(name: str) -> None:
    if DESCRIPTIONS.spec(name).stata is None:
        pytest.skip("an addition of the package: no Stata command")
    cases = "\n".join(
        path.read_text(encoding="utf-8") for path in (ROOT / "validation" / "cases").glob("*.py")
    )
    assert name in cases, f"{name!r} claims a Stata command but no validation case runs it"


@pytest.mark.parametrize("name", KINDS.names())
def test_every_figure_has_a_gallery_picture_and_is_documented(name: str) -> None:
    assert (ROOT / "docs" / "gallery" / f"{name}.png").exists(), f"no picture for {name!r}"
    documentation = (ROOT / "docs" / "documentation.md").read_text(encoding="utf-8")
    assert f'`"{name}"`' in documentation, f"figure {name!r} is not named in part 6"


@pytest.mark.parametrize("name", CHECKS.names())
def test_every_alert_code_of_a_check_has_a_message_and_an_advice_in_both_languages(
    name: str,
) -> None:
    for code in CHECKS.spec(name).codes:
        for key in (f"explore.audit.{code}", f"explore.audit.{code}.advice"):
            for lang in ("en", "fr"):
                assert _i18n.has_language(key, lang), f"{key}: {lang} text missing"
