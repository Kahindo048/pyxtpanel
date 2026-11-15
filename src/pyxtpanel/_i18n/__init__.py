"""Message catalogue (French and English), one table per domain.

Layer 0.  Every module that needs messages declares them next to its code with
:func:`register`, under a domain prefix (``panel.duplicates``, ``diag.ips.too_few_periods``),
so adding a feature never edits a shared file.  A contract test imports every module and
checks that each key exists in both languages with the same fields.
"""

from __future__ import annotations

from collections.abc import Mapping
from string import Formatter

from pyxtpanel._config import get_options
from pyxtpanel._types import Lang

_CATALOGUE: dict[str, dict[str, str]] = {}


def register(domain: str, table: Mapping[str, Mapping[str, str]]) -> None:
    """Add messages to the catalogue.

    Parameters
    ----------
    domain
        Prefix every key of ``table`` must start with, followed by a dot.
    table
        ``{key: {"en": text, "fr": text}}``; texts are ``str.format`` templates.

    Raises
    ------
    ValueError
        If a key does not start with ``domain + "."`` or is already registered with a
        different text (registering the same table twice is harmless).
    """
    for key, texts in table.items():
        if not key.startswith(domain + "."):
            raise ValueError(f"message key {key!r} must start with {domain + '.'!r}")
        if key in _CATALOGUE and _CATALOGUE[key] != dict(texts):
            raise ValueError(f"message key {key!r} is already registered")
        _CATALOGUE[key] = dict(texts)


def message(key: str, lang: Lang | None = None, **fields: object) -> str:
    """Return the message registered under ``key`` in the requested language.

    Parameters
    ----------
    key
        Message identifier.
    lang
        ``"en"`` or ``"fr"``; defaults to the global ``lang`` option.  A key missing in the
        requested language falls back to English.
    **fields
        Values substituted in the template.

    Raises
    ------
    KeyError
        If the key is not registered.
    """
    language: Lang = lang if lang is not None else get_options().lang
    texts = _CATALOGUE[key]
    return (texts.get(language) or texts["en"]).format(**fields)


def all_keys() -> list[str]:
    """Return all registered message keys, sorted."""
    return sorted(_CATALOGUE)


def has_language(key: str, lang: Lang) -> bool:
    """Whether ``key`` is translated into ``lang``."""
    return lang in _CATALOGUE[key]


def fields_of(key: str, lang: Lang) -> frozenset[str]:
    """Names of the fields of the template of ``key`` in ``lang``."""
    return frozenset(
        name for _, name, _, _ in Formatter().parse(_CATALOGUE[key][lang]) if name is not None
    )
