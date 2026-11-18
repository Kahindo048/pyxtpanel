"""Metadata and registry of the statistical tests (L.1.7).

A test is a pure function that returns a :class:`~pyxtpanel.diagnostics._core.result.TestResult`
and is registered with :func:`register_test`.  Its :class:`TestSpec` is the single source for the
list of tests, the documentation table, the completeness checks of ``tests/contract/`` and the
validation coverage table.

**Validation against Stata is not optional**: ``TestSpec.validation`` lists the cases of
``validation/cases/diagnostics_<family>.py`` that were run in Stata and frozen in
``tests/golden/``, and a contract test fails for a test without one.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, TypeVar

import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._registry import Registry

if TYPE_CHECKING:
    from pyxtpanel.diagnostics._core.result import TestResult

F = TypeVar("F", bound=Callable[..., object])

#: Families of tests, one folder of ``diagnostics/`` each.
FAMILIES = (
    "unitroot",
    "crosssection",
    "cointegration",
    "serial",
    "causality",
    "postestimation",
    "heteroskedasticity",
    "effects",
    "normality",
    "specification",
    "heterogeneity",
)


@dataclass(frozen=True)
class TestSpec:
    """What a test is, what it needs and how it is validated.

    Attributes
    ----------
    family
        Folder of ``diagnostics/`` the test lives in (one of :data:`FAMILIES`).
    category
        ``"A"`` for a test on the data (input: a panel and variable names) or ``"B"`` for a
        generic post-estimation test (input: bare arrays).  Tests that need a model object are
        category C and do not belong here.
    wave
        Delivery wave of the plan (1 most used, 2 complements, 3 rarer tests).
    reference
        Publication the statistic comes from.
    stata
        The Stata command the test is validated against (``"xtunitroot ips"``); ``"mata"`` when
        Stata has no command and the case compares with an independent Mata implementation of
        the published formula.  Never empty.
    validation
        Tags of the validation cases (``validation/cases/diagnostics_<family>.py``) run in Stata
        and frozen in ``tests/golden/``.  Never empty.
    inputs
        ``"panel"`` for category A, ``"arrays"`` for category B.
    balanced
        ``"required"`` when the test needs a strongly balanced panel, ``"both"`` when it also
        accepts an unbalanced one, ``"na"`` for arrays.
    r
        The equivalent in R (``"plm::purtest(test='ips')"``), ``None`` when there is none.
    tables
        Names of the embedded tables the test reads (``_data/tables/``).
    aliases
        Other names the test answers to, in lower case: the names of the Stata subcommands
        (``"ips"``).  An alias is never a test name nor the alias of another test.
    summary
        One-line description.
    """

    # Not a pytest test class despite the name.
    __test__ = False

    family: str
    category: Literal["A", "B"]
    wave: Literal[1, 2, 3]
    reference: str
    stata: str
    validation: tuple[str, ...]
    inputs: Literal["panel", "arrays"]
    summary: str
    balanced: Literal["required", "both", "na"] = "na"
    r: str | None = None
    tables: tuple[str, ...] = ()
    aliases: tuple[str, ...] = ()


_i18n.register(
    "diag",
    {
        "diag.dispatch.unknown_test": {
            "en": "unknown test {name!r}; available: {available}.",
            "fr": "test inconnu {name!r} ; disponibles : {available}.",
        },
    },
)


#: Registry of the tests (``@register_test("ips", TestSpec(...))``).
TESTS: Registry[TestSpec] = Registry("test", package="pyxtpanel.diagnostics")


def register_test(name: str, spec: TestSpec) -> Callable[[F], F]:
    """Register a statistical test under ``name``.

    Parameters
    ----------
    name
        Registered name of the test.
    spec
        Its :class:`TestSpec`: family, category, inputs and references.

    Returns
    -------
    callable
        The decorator that registers the function.
    """
    return TESTS.register(name, spec)


def list_tests() -> pd.DataFrame:
    """Return the registered tests, one row each, with their metadata.

    Columns: ``family``, ``category``, ``wave``, ``inputs``, ``balanced``, ``stata``, ``r``,
    ``reference``, ``summary``; indexed by test name.

    Returns
    -------
    pandas.DataFrame
        One row per test, indexed by name.
    """
    rows = {
        name: {
            "family": spec.family,
            "category": spec.category,
            "wave": spec.wave,
            "inputs": spec.inputs,
            "balanced": spec.balanced,
            "stata": spec.stata,
            "r": spec.r,
            "reference": spec.reference,
            "summary": spec.summary,
        }
        for name, spec in TESTS.specs()
    }
    columns = ["family", "category", "wave", "inputs", "balanced", "stata", "r", "reference"]
    frame = pd.DataFrame.from_dict(rows, orient="index", columns=[*columns, "summary"])
    return frame.rename_axis("test")


def resolve_test(name: str) -> str:
    """Return the registered name of the test called ``name`` (a name or an alias).

    Parameters
    ----------
    name
        Name or alias of the test.

    Returns
    -------
    str
        The registered name.

    Raises
    ------
    SpecificationError
        If no test has that name or alias; the message lists the available ones.
    """
    key = name.lower()
    names = TESTS.names()
    if key in names:
        return key
    for registered, spec in TESTS.specs():
        if key in spec.aliases:
            return registered
    aliases = sorted(alias for _, spec in TESTS.specs() for alias in spec.aliases)
    raise SpecificationError(
        _i18n.message(
            "diag.dispatch.unknown_test", name=name, available=", ".join([*names, *aliases])
        ),
        key="diag.dispatch.unknown_test",
    )


def run_test(name: str, *args: object, **kwargs: object) -> TestResult:
    """Run the registered test ``name`` (or one of its aliases) with the given arguments.

    ``run_test("ips", panel, "x", lags=1)`` is ``im_pesaran_shin(panel, "x", lags=1)``.

    Parameters
    ----------
    name
        Registered name of the test, or one of its aliases.
    args
        Arguments of the test (the panel or the arrays, then the variable names).
    kwargs
        Options of the test.

    Returns
    -------
    TestResult
        The test statistic, its p-value and the quantities behind them.

    Raises
    ------
    SpecificationError
        If no test has that name; the message lists the available names.
    """
    result: TestResult = TESTS.get(resolve_test(name))(*args, **kwargs)  # type: ignore[assignment]
    return result
