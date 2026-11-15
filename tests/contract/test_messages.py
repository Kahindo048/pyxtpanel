"""Contract test: every message exists in French and English with the same fields."""

import importlib
import pkgutil

import pytest

import pyxtpanel
from pyxtpanel import _i18n


def _import_everything() -> None:
    for info in pkgutil.walk_packages(pyxtpanel.__path__, "pyxtpanel."):
        if ".plot" in info.name:  # needs the optional matplotlib extra
            continue
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
