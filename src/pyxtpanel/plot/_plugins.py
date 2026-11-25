"""Figure kinds supplied by other packages: the entry-point group ``pyxtpanel.plot_kinds``.

A companion package (``pyxtviz``) declares, in its ``pyproject.toml``::

    [project.entry-points."pyxtpanel.plot_kinds"]
    pyxtviz = "pyxtviz.plugin:KINDS"

The target is an object with two methods: ``names()``, the kinds it supplies, and ``get(name)``,
the function that draws one (``f(panel, ..., lang=None) -> Figure``).  Rules:

* **Lazy.**  Nothing is looked up unless the kind asked for is not one of the package's own, so
  ``import pyxtpanel`` and the built-in figures never pay for it.
* **The package's own kinds win** when a plugin supplies the same name.
* Among plugins, the one whose entry point is named first (alphabetical order) wins.
* **A plugin that cannot be loaded** never breaks the built-in figures; it raises a clear error
  only when the kind asked for was not found elsewhere (it may have been in that plugin).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any, Protocol

from pyxtpanel import _i18n
from pyxtpanel._errors import PyxtpanelError
from pyxtpanel._types import Lang

#: Name of the entry-point group.
GROUP = "pyxtpanel.plot_kinds"

_i18n.register(
    "plot",
    {
        "plot.plugin_broken": {
            "en": "the figure plugin {plugin!r} ({target}) cannot be used: {reason}. Reinstall or "
            "upgrade it, or uninstall it.",
            "fr": "le greffon de figures {plugin!r} ({target}) est inutilisable : {reason}. "
            "Réinstallez-le ou mettez-le à jour, ou désinstallez-le.",
        },
        "plot.unknown_kind": {
            "en": "unknown figure kind {kind!r}; available: {available}.",
            "fr": "type de figure inconnu {kind!r} ; disponibles : {available}.",
        },
        "plot.unknown_kind_install": {
            "en": "unknown figure kind {kind!r}; available: {available}. More figures come with "
            "the companion package pyxtviz: pip install pyxtviz.",
            "fr": "type de figure inconnu {kind!r} ; disponibles : {available}. D'autres figures "
            "viennent avec le paquet compagnon pyxtviz : pip install pyxtviz.",
        },
    },
)


class KindSource(Protocol):
    """What an entry point of the group must point to."""

    def names(self) -> Iterable[str]:
        """Return the names of the kinds supplied."""
        ...

    def get(self, name: str) -> Callable[..., Any]:
        """Return the function that draws the kind ``name``."""
        ...


_cache: list[Any] | None = None


def _discover() -> list[Any]:
    """Entry points of the group, sorted by name (the scan of the installed packages)."""
    from importlib.metadata import entry_points

    return sorted(entry_points(group=GROUP), key=lambda entry: entry.name)


def _entries() -> list[Any]:
    global _cache
    if _cache is None:
        _cache = _discover()
    return _cache


def _reset() -> None:
    """Forget what was discovered (for tests that install or remove a plugin)."""
    global _cache
    _cache = None


def _broken(entry: Any, error: BaseException, lang: Lang | None) -> PyxtpanelError:
    return PyxtpanelError(
        _i18n.message(
            "plot.plugin_broken",
            lang,
            plugin=entry.name,
            target=getattr(entry, "value", "?"),
            reason=f"{type(error).__name__}: {error}",
        ),
        key="plot.plugin_broken",
    )


def has_plugins() -> bool:
    """Whether some package supplies figure kinds (loaded or not)."""
    return bool(_entries())


def _listed(entry: Any) -> set[str] | None:
    """Lower-case names a plugin supplies, or ``None`` if it cannot be loaded."""
    try:
        return {str(name).lower() for name in entry.load().names()}
    except Exception:  # noqa: BLE001 - any failure of a third-party plugin
        return None


def plugin_names() -> list[str]:
    """Return the names of the kinds the working plugins supply, sorted (broken ones skipped)."""
    found: set[str] = set()
    for entry in _entries():
        found |= _listed(entry) or set()
    return sorted(found)


def find_kind(kind: str, lang: Lang | None = None) -> Callable[..., Any] | None:
    """Return the function drawing ``kind`` from the first plugin that supplies it.

    Parameters
    ----------
    kind
        Name of the figure kind (case does not matter).
    lang
        ``"en"`` or ``"fr"``: language of the error, if any (default: the global option).

    Returns
    -------
    callable or None
        The function, or ``None`` when no plugin supplies the kind.

    Raises
    ------
    PyxtpanelError
        If no working plugin supplies the kind and some plugin could not be loaded, or if the
        plugin that lists the kind fails to give it.
    """
    wanted = str(kind).lower()
    failure: tuple[Any, BaseException] | None = None
    for entry in _entries():
        try:
            source = entry.load()
            listed = {str(name).lower() for name in source.names()}
        except Exception as error:  # noqa: BLE001 - any failure of a third-party plugin
            failure = failure or (entry, error)
            continue
        if wanted not in listed:
            continue
        try:
            draw = source.get(wanted)
        except Exception as error:  # noqa: BLE001
            raise _broken(entry, error, lang) from error
        if not callable(draw):
            raise _broken(entry, TypeError("get() did not return a function"), lang)
        return draw
    if failure is not None:
        raise _broken(failure[0], failure[1], lang) from failure[1]
    return None
