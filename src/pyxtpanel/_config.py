"""Global options of the package (layer 0).

This is the only module holding global state.  Options are read when a function runs, never
at import time.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, fields, replace
from typing import Any

from pyxtpanel._types import LANGUAGES, Lang


@dataclass(frozen=True)
class Options:
    """Immutable set of global options.

    Attributes
    ----------
    lang : {"en", "fr"}
        Language of messages and reports.
    """

    lang: Lang = "en"


_options = Options()


def get_options() -> Options:
    """Return the current global options."""
    return _options


def set_options(**changes: Any) -> Options:
    """Change global options and return the previous ones.

    Parameters
    ----------
    **changes
        Any field of :class:`Options`: ``lang`` (``"en"`` or ``"fr"``).

    Returns
    -------
    Options
        The options before the change, to restore them later.

    Raises
    ------
    ValueError
        If an option name or value is invalid.

    Examples
    --------
    >>> previous = set_options(lang="fr")
    >>> previous.lang
    'en'
    >>> _ = set_options(lang=previous.lang)
    """
    global _options
    previous = _options
    unknown = set(changes) - {f.name for f in fields(Options)}
    if unknown:
        raise ValueError(f"Unknown option(s): {sorted(unknown)}")
    candidate = replace(previous, **changes)
    if candidate.lang not in LANGUAGES:
        raise ValueError(f"lang must be one of {list(LANGUAGES)}, got {candidate.lang!r}")
    _options = candidate
    return previous


@contextmanager
def option_context(**changes: Any) -> Iterator[Options]:
    """Change options inside a ``with`` block and restore them afterwards.

    Examples
    --------
    >>> with option_context(lang="fr") as options:
    ...     options.lang
    'fr'
    >>> get_options().lang
    'en'
    """
    previous = set_options(**changes)
    try:
        yield get_options()
    finally:
        set_options(**{f.name: getattr(previous, f.name) for f in fields(Options)})
