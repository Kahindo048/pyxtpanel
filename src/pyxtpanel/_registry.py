"""Generic registry behind every extension point of the package (layer 0).

A feature (a statistical test, an operator, a descriptive statistic, an audit check, a plot
kind, an output format) is one module that registers itself with a decorator.  The registry
keeps its metadata (a frozen dataclass specific to the extension point), so that lists, the
documentation and the completeness checks are built from one source.
"""

from __future__ import annotations

import importlib
import pkgutil
from collections.abc import Callable, Iterator
from typing import Generic, TypeVar

from pyxtpanel._errors import SpecificationError

S = TypeVar("S")
F = TypeVar("F", bound=Callable[..., object])


class Registry(Generic[S]):
    """Named entries, each a callable with its metadata.

    Parameters
    ----------
    kind
        Name of what is registered (``"test"``, ``"operator"``...), used in error messages.
    package
        Dotted name of the package whose submodules hold the entries.  :meth:`load_all`
        imports them so that their decorators run.
    """

    def __init__(self, kind: str, package: str | None = None) -> None:
        self.kind = kind
        self.package = package
        self._entries: dict[str, tuple[S, Callable[..., object]]] = {}
        self._loaded = False

    def register(self, name: str, spec: S) -> Callable[[F], F]:
        """Return a decorator registering a function under ``name`` with metadata ``spec``.

        Raises
        ------
        SpecificationError
            If ``name`` is already registered by another function.
        """
        key = name.lower()

        def decorator(func: F) -> F:
            known = self._entries.get(key)
            if known is not None and known[1].__qualname__ != func.__qualname__:
                raise SpecificationError(
                    f"{self.kind} {key!r} is already registered", key="registry.duplicate"
                )
            self._entries[key] = (spec, func)
            return func

        return decorator

    def load_all(self) -> None:
        """Import every submodule of :attr:`package` (once), running their decorators."""
        if self._loaded or self.package is None:
            return
        self._loaded = True
        root = importlib.import_module(self.package)
        for info in pkgutil.walk_packages(root.__path__, root.__name__ + "."):
            importlib.import_module(info.name)

    def get(self, name: str) -> Callable[..., object]:
        """Return the function registered under ``name``.

        Raises
        ------
        SpecificationError
            If no such entry exists; the message lists the available names.
        """
        self.load_all()
        try:
            return self._entries[name.lower()][1]
        except KeyError:
            raise SpecificationError(
                f"unknown {self.kind} {name!r}; available: {', '.join(self.names())}",
                key="registry.unknown",
            ) from None

    def spec(self, name: str) -> S:
        """Return the metadata registered under ``name``."""
        self.get(name)
        return self._entries[name.lower()][0]

    def names(self) -> list[str]:
        """Return the registered names, sorted."""
        self.load_all()
        return sorted(self._entries)

    def specs(self) -> Iterator[tuple[str, S]]:
        """Iterate over ``(name, spec)`` pairs in name order."""
        for name in self.names():
            yield name, self._entries[name][0]

    def __contains__(self, name: object) -> bool:
        self.load_all()
        return isinstance(name, str) and name.lower() in self._entries

    def __len__(self) -> int:
        self.load_all()
        return len(self._entries)
