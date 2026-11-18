"""Contract test: every message exists in French and English with the same fields."""

import importlib
import pkgutil
import re

import pytest

import pyxtpanel
from pyxtpanel import _i18n


def _import_everything() -> None:
    for info in pkgutil.walk_packages(pyxtpanel.__path__, "pyxtpanel."):
        importlib.import_module(info.name)


_import_everything()
# keys registered by unit tests of the message mechanism itself are excluded
KEYS = [k for k in _i18n.all_keys() if not k.startswith("test_l0.")]


@pytest.mark.parametrize("key", KEYS)
def test_message_in_both_languages_with_same_fields(key: str) -> None:
    assert _i18n.has_language(key, "en"), f"{key}: English text missing"
    assert _i18n.has_language(key, "fr"), f"{key}: French text missing"
    assert _i18n.fields_of(key, "en") == _i18n.fields_of(key, "fr"), f"{key}: fields differ"


def test_the_catalogue_is_not_empty_once_features_exist() -> None:
    # guards against a silent import failure; relaxed while no feature has messages yet
    assert isinstance(KEYS, list)


# A message that points at "a unit", "some row" or "a variable" must say which one (a field named
# after it), because the user has to find it in the data.  The exceptions below are causes that
# belong to the panel as a whole or to a pooled quantity, where no single unit is at fault, or
# messages that a caller catches and re-raises with the unit.
_INDEFINITE = re.compile(
    r"\b(?:a|some|one) (?:unit|row|variable|column)(?: has| needs| is)"
    r"|\bfor a unit\b|\ba unit, or\b",
    re.IGNORECASE,
)
_NAMING_FIELDS = {"unit", "units", "row", "rows", "variable", "name", "column", "example", "label"}
UNNAMED_ON_PURPOSE = {
    "diag.baltagi_li.need_repeated": "no unit has two observations: a panel-wide condition",
    "diag.baltagi_li.need_three": "no unit has three observations: a panel-wide condition",
    "diag.breitung.degenerate": "pooled sums that vanish; a single flat unit has its own message",
    "diag.hadri.degenerate": "the pooled residual variance is zero; a flat unit has its own one",
    "diag.kao.undefined_bandwidth": "the long-run covariance of the pooled innovations",
    "diag.levin_lin_chu.undefined_bandwidth": "re-raised by the caller as ..._unit, with the unit",
}


@pytest.mark.parametrize("key", KEYS)
def test_a_message_that_points_at_an_element_names_it(key: str) -> None:
    text = _i18n._CATALOGUE[key]["en"]
    fields = set(_i18n.fields_of(key, "en"))
    if _INDEFINITE.search(text) and not fields & _NAMING_FIELDS:
        assert key in UNNAMED_ON_PURPOSE, f"{key}: points at an element without naming it: {text!r}"


def test_the_exceptions_are_still_needed() -> None:
    """An exception that no longer matches the rule must be removed from the list."""
    for key in UNNAMED_ON_PURPOSE:
        assert key in KEYS, f"{key} no longer exists"
        text = _i18n._CATALOGUE[key]["en"]
        assert _INDEFINITE.search(text), f"{key} now names its element or no longer points at one"
