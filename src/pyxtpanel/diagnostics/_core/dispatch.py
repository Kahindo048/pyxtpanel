"""Run a category A test on a declared panel by name (L.1.6).

``unit_root(panel, "y", test="ips", lags=1)``, ``PanelData.test("ips", "y", lags=1)`` and the
estimation packages all end here.  The name is resolved (names and Stata aliases), the test is
checked to be a test on the data and, when a family is given, to belong to it, and the call is
checked against the signature of the test before it runs: an unknown option or a missing
argument is reported in French or English with the options the test accepts, instead of Python's
bare ``TypeError``.  Nothing is ever ignored silently, and no test is chosen by default: the
tests of a family answer different questions.
"""

from __future__ import annotations

import inspect

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.diagnostics._core.result import TestResult
from pyxtpanel.diagnostics._core.spec import TESTS, resolve_test
from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.dispatch.test_required": {
            "en": "choose a test with test=...; for {scope}: {available}.",
            "fr": "choisissez un test avec test=... ; pour {scope} : {available}.",
        },
        "diag.dispatch.not_panel_test": {
            "en": "{name!r} is a test on arrays (category B), not on a panel: call "
            "pyxtpanel.diagnostics.{name}(...) or run_test({name!r}, ...) with the arrays.",
            "fr": "{name!r} est un test sur tableaux (catégorie B), pas sur un panel : appelez "
            "pyxtpanel.diagnostics.{name}(...) ou run_test({name!r}, ...) avec les tableaux.",
        },
        "diag.dispatch.wrong_family": {
            "en": "{name!r} is not a {family} test; the {family} tests are: {available}.",
            "fr": "{name!r} n'est pas un test de la famille {family} ; ceux de cette famille "
            "sont : {available}.",
        },
        "diag.dispatch.bad_call": {
            "en": "{name}: {reason}. Arguments of the test: {accepted}.",
            "fr": "{name} : {reason}. Arguments du test : {accepted}.",
        },
        "diag.dispatch.not_declared": {
            "en": "a test on the data needs a declared panel (use xtset(frame, id, time)), got "
            "{kind}.",
            "fr": "un test sur les données exige un panel déclaré (utilisez xtset(frame, id, "
            "temps)), reçu {kind}.",
        },
    },
)

_FAMILY_LABELS = {
    "unitroot": "unit-root",
    "crosssection": "cross-section dependence",
    "serial": "serial-correlation",
    "heterogeneity": "slope-heterogeneity",
}


def panel_tests(family: str | None = None) -> list[str]:
    """Return the names of the tests on a panel (category A), optionally of one family."""
    return [
        name
        for name, spec in TESTS.specs()
        if spec.category == "A" and (family is None or spec.family == family)
    ]


def run_panel_test(
    panel: object,
    test: str | None,
    *args: object,
    family: str | None = None,
    **options: object,
) -> TestResult:
    """Run the category A test ``test`` on ``panel`` with the given arguments and options.

    Parameters
    ----------
    panel
        A declared panel.
    test
        Name or alias of the test (``"ips"``, ``"im_pesaran_shin"``); required.
    *args
        Positional arguments after the panel (the variable names).
    family
        Restrict the test to a family (``"unitroot"``, ``"crosssection"``, ``"serial"``,
        ``"heterogeneity"``).
    **options
        The options of the test.

    Raises
    ------
    SpecificationError
        If the test is missing, unknown, not a test on a panel, outside the family, or called
        with an argument or an option it does not accept.
    """
    scope = _FAMILY_LABELS.get(family or "", "panel")
    candidates = ", ".join(panel_tests(family))
    if test is None:
        raise SpecificationError(
            _i18n.message("diag.dispatch.test_required", scope=scope, available=candidates),
            key="diag.dispatch.test_required",
        )
    name = resolve_test(test)
    spec = TESTS.spec(name)
    if spec.category != "A":
        raise SpecificationError(
            _i18n.message("diag.dispatch.not_panel_test", name=name),
            key="diag.dispatch.not_panel_test",
        )
    if family is not None and spec.family != family:
        raise SpecificationError(
            _i18n.message(
                "diag.dispatch.wrong_family", name=name, family=scope, available=candidates
            ),
            key="diag.dispatch.wrong_family",
        )
    if not isinstance(panel, PanelBase):
        raise SpecificationError(
            _i18n.message("diag.dispatch.not_declared", kind=type(panel).__name__),
            key="diag.dispatch.not_declared",
        )
    function = TESTS.get(name)
    signature = inspect.signature(function)
    try:
        signature.bind(panel, *args, **options)
    except TypeError as error:
        raise SpecificationError(
            _i18n.message(
                "diag.dispatch.bad_call",
                name=name,
                reason=str(error),
                accepted=_accepted(signature),
            ),
            key="diag.dispatch.bad_call",
        ) from None
    result = function(panel, *args, **options)
    assert isinstance(result, TestResult)
    return result


def _accepted(signature: inspect.Signature) -> str:
    """Return the arguments after the panel, keyword-only ones marked with their defaults."""
    parts = []
    for index, (name, parameter) in enumerate(signature.parameters.items()):
        if index == 0:
            continue  # the panel
        if parameter.default is inspect.Parameter.empty:
            parts.append(name)
        else:
            parts.append(f"{name}={parameter.default!r}")
    return ", ".join(parts)
