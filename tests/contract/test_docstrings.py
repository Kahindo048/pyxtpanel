"""The docstrings of the public API are complete (J.2.6, M.3).

``tools/docstrings.py`` lists the public objects (the names exported by the package and its public
sub-packages, the methods of ``PanelData``, every registered test and figure kind) and checks the
content of each docstring: a summary, every parameter described, a ``Returns`` (or ``Yields``)
section, a ``Raises`` section when the function raises an error of the package.  The first test is
the gate; the others check that the checker finds what it is meant to find.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest
from tools import docstrings

from pyxtpanel._errors import SpecificationError


def test_the_public_api_has_no_docstring_defect() -> None:
    count, problems = docstrings.check_all()
    assert count > 150, "the public objects were not found"
    assert problems == [], "\n".join(problems)


def complete(x: int, *, y: int = 0) -> int:
    """Add two numbers.

    Parameters
    ----------
    x
        The first.
    y
        The second.

    Returns
    -------
    int
        The sum.
    """
    return x + y


def test_a_complete_docstring_passes() -> None:
    assert docstrings.check_function("complete", complete) == []


def test_a_parameter_without_description_is_found() -> None:
    def f(x: int, y: int) -> None:
        """Do something.

        Parameters
        ----------
        x
            The first.
        """

    assert docstrings.check_function("f", f) == ["f: parameter 'y' is not documented"]


def test_a_stale_parameter_is_found() -> None:
    def f(x: int) -> None:
        """Do something.

        Parameters
        ----------
        x, old
            The first, and a parameter that was renamed.
        """

    assert docstrings.check_function("f", f) == ["f: documented parameter 'old' does not exist"]


def test_the_returns_section_is_required_when_something_is_returned() -> None:
    def f() -> int:
        """Return a number."""
        return 1

    assert docstrings.check_function("f", f) == ["f: no Returns section"]


def test_a_generator_needs_yields_not_returns() -> None:
    def g() -> Iterator[int]:
        """Count.

        Returns
        -------
        int
            Wrong section for a generator.
        """
        yield 1

    assert docstrings.check_function("g", g) == ["g: no Yields section"]


def test_an_error_of_the_package_must_be_documented() -> None:
    def f(x: int) -> None:
        """Check ``x``.

        Parameters
        ----------
        x
            The value.
        """
        if x < 0:
            raise SpecificationError("negative")

    assert docstrings.check_function("f", f) == [
        "f: raises an error of the package but has no Raises section"
    ]


def test_a_missing_docstring_is_found() -> None:
    def f() -> None:
        pass

    assert docstrings.check_function("f", f) == ["f: no docstring"]


@pytest.mark.parametrize(
    ("body", "names"),
    [
        ("a\n    One.\nb, c : int\n    Two.\n**options\n    Three.", {"a", "b", "c", "options"}),
        ("", set()),
    ],
)
def test_the_parameter_section_is_read(body: str, names: set[str]) -> None:
    assert docstrings.documented_parameters(body) == names
