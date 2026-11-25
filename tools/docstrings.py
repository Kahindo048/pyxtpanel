"""Check the docstrings of the public API (J.2.6, M.3).

    python -m tools.docstrings           # list the defects, exit code 1 if there is any
    python -m tools.docstrings --count   # only the counts

What is checked, for every name of the public API (``__all__`` of the package and of its public
sub-packages, the methods of ``PanelData``, every registered test and figure kind):

* a docstring with a one-line summary;
* every parameter of a function or method is described in the ``Parameters`` section, and the
  section describes nothing that is not a parameter (a renamed argument leaves a stale entry);
* a function that returns something has a ``Returns`` section (``Yields`` for a generator);
* a function that raises an exception of the package itself (``SpecificationError``,
  ``DataError``...) has a ``Raises`` section.

The format is numpydoc, the one the documentation uses; ``ruff`` checks the layout (summary line,
blank lines, section underlines), this script checks the content.
"""

from __future__ import annotations

import ast
import importlib
import inspect
import re
import sys
import textwrap
from collections.abc import Callable, Iterator
from typing import Any

#: Modules whose ``__all__`` is public API.
MODULES = (
    "pyxtpanel",
    "pyxtpanel.contract",
    "pyxtpanel.diagnostics",
    "pyxtpanel.extend",
    "pyxtpanel.plot",
)

_ERRORS = {
    "SpecificationError",
    "DataError",
    "DeclarationError",
    "TimeAxisError",
    "NumericalError",
    "PyxtpanelError",
}
_HEADER = re.compile(r"^(\w[\w ]*)\n-{3,}\s*$", re.MULTILINE)


def sections(doc: str) -> dict[str, str]:
    """Split a numpydoc docstring into ``{section title: body}`` (the summary is ``""``)."""
    text = inspect.cleandoc(doc)
    found = list(_HEADER.finditer(text))
    parts = {"": text[: found[0].start()] if found else text}
    for index, match in enumerate(found):
        end = found[index + 1].start() if index + 1 < len(found) else len(text)
        parts[match.group(1).strip()] = text[match.end() : end]
    return parts


def documented_parameters(body: str) -> set[str]:
    """Return the names described in a ``Parameters`` section (``a, b : type`` gives both)."""
    names: set[str] = set()
    for line in body.splitlines():
        if not line.strip() or line[0] in " \t":
            continue  # blank or description (indented)
        for name in line.split(":", 1)[0].split(","):
            names.add(name.strip().lstrip("*"))
    return names


def _parameters(function: Callable[..., Any]) -> list[str]:
    return [name for name in inspect.signature(function).parameters if name not in ("self", "cls")]


def _raises_a_package_error(function: Callable[..., Any]) -> bool:
    try:
        tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
    except (OSError, TypeError, SyntaxError):
        return False
    for node in ast.walk(tree):
        if isinstance(node, ast.Raise) and node.exc is not None:
            target = node.exc.func if isinstance(node.exc, ast.Call) else node.exc
            if isinstance(target, ast.Name) and target.id in _ERRORS:
                return True
    return False


def _returns_something(function: Callable[..., Any]) -> bool:
    annotation = inspect.signature(function).return_annotation
    return annotation not in (inspect.Signature.empty, None, "None", "NoReturn", "Never")


def check_function(label: str, function: Callable[..., Any]) -> list[str]:
    """Return the defects of the docstring of a function or method."""
    doc = inspect.getdoc(function)
    if not doc or not doc.strip():
        return [f"{label}: no docstring"]
    parts = sections(doc)
    problems = []
    summary = parts[""].strip().splitlines()
    if not summary or not summary[0].strip():
        problems.append(f"{label}: empty summary")
    documented = documented_parameters(parts.get("Parameters", ""))
    expected = set(_parameters(function))
    for name in sorted(expected - documented):
        problems.append(f"{label}: parameter {name!r} is not documented")
    for name in sorted(documented - expected):
        problems.append(f"{label}: documented parameter {name!r} does not exist")
    wants_yields = inspect.isgeneratorfunction(inspect.unwrap(function))
    if wants_yields and "Yields" not in parts:
        problems.append(f"{label}: no Yields section")
    elif not wants_yields and _returns_something(function) and "Returns" not in parts:
        problems.append(f"{label}: no Returns section")
    if _raises_a_package_error(function) and "Raises" not in parts:
        problems.append(f"{label}: raises an error of the package but has no Raises section")
    return problems


def check_class(label: str, cls: type) -> list[str]:
    """Return the defects of the docstring of a class (a summary is enough)."""
    doc = inspect.getdoc(cls)
    if not doc or not doc.strip().splitlines()[0].strip():
        return [f"{label}: no docstring"]
    return []


def _exported() -> Iterator[tuple[str, Any]]:
    """Yield the functions and classes the public modules export."""
    for module_name in MODULES:
        module = importlib.import_module(module_name)
        for name in module.__all__:
            value = getattr(module, name)
            if inspect.isfunction(value) or inspect.isclass(value):
                yield f"{module_name}.{name}", value


def _registered() -> Iterator[tuple[str, Any]]:
    """Yield the methods of ``PanelData``, the registered tests and the figure kinds."""
    from pyxtpanel import PanelData
    from pyxtpanel.diagnostics._core.spec import TESTS
    from pyxtpanel.plot import KINDS

    for name, method in inspect.getmembers(PanelData, inspect.isfunction):
        if not name.startswith("_"):
            yield f"PanelData.{name}", method
    for name in TESTS.names():
        yield f"test {name!r}", TESTS.get(name)
    for name in KINDS.names():
        yield f"figure {name!r}", KINDS.get(name)


def public_objects() -> Iterator[tuple[str, Any]]:
    """Yield ``(label, object)`` for each documented entry point of the package, once."""
    seen: set[int] = set()
    for label, obj in (*_exported(), *_registered()):
        if id(obj) not in seen:
            seen.add(id(obj))
            yield label, obj


def check_all() -> tuple[int, list[str]]:
    """Return the number of objects checked and the defects."""
    problems: list[str] = []
    count = 0
    for label, obj in public_objects():
        count += 1
        problems += check_class(label, obj) if inspect.isclass(obj) else check_function(label, obj)
    return count, problems


def main(argv: list[str]) -> int:
    """Print the defects; the exit code is 1 when there is any."""
    count, problems = check_all()
    if "--count" not in argv:
        for problem in problems:
            print(problem)
    print(f"{count} public objects checked, {len(problems)} defects")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
