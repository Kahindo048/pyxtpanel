"""Create the files of a new feature from a template (A.5.9).

    python -m tools.nouveau operateur  moving_median
    python -m tools.nouveau description xtgini
    python -m tools.nouveau controle    short_units
    python -m tools.nouveau graphique   heatmap
    python -m tools.nouveau test        fisher --famille unitroot --categorie A

Creates, without touching any existing file:

* the module (``ops/<name>.py``, ``describe/<name>.py``, ``explore/checks/<name>.py`` or
  ``plot/kinds/<name>.py``) that registers itself,
* its unit test (``tests/unit/test_<package>_<name>.py``, ``tests/plot/test_figure_<name>.py``
  for a figure),
* a skeleton of validation case (``validation/cases/<name>.py``) when the feature has a Stata
  reference.

The registry finds the module by itself (``Registry.load_all``), so nothing else has to be
edited, except two lines the tool reminds you of: the sample call of a description in
``tests/contract/test_public_api.py`` and the entry of ``CHANGELOG.md``.

A statistical test (``test``) gets its module in ``diagnostics/<family>/``, its unit test in
``tests/unit/diagnostics/`` and its Stata validation case in
``validation/cases/diagnostics_<family>_<name>.py`` (domain ``diagnostics_<family>``).  The
contract test of ``tests/contract/`` fails until the case has been run in Stata and frozen.

A figure (``graphique``) gets ``plot/kinds/<name>.py``, which holds a pure function that computes
what is drawn (``<name>_data``, tested without matplotlib) and the drawing function registered as
a kind, and ``tests/plot/test_figure_<name>.py``.  The figure must then be named in part 6 of
``docs/documentation.md`` (a contract test fails otherwise) and gets its gallery picture with
``python -m tools.galerie``.

There is no template for an output format: the formats of a result (``to_markdown``, ``to_latex``,
HTML) are methods of :class:`~pyxtpanel._results.Result` written with the helpers of
``pyxtpanel/_format.py``, and no registry holds them, so a new format is one helper there and one
method (decision of 2026-11-18, A.5.9).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KINDS = ("operateur", "description", "controle", "graphique", "test")

OPERATOR = '''"""{title} (operator or transformation)."""

from __future__ import annotations

from pyxtpanel._types import FloatArray
from pyxtpanel.ops._spec import OpSpec, operator
from pyxtpanel.panel.index import PanelIndex


@operator(
    "{name}",
    OpSpec(
        kind="transform",  # "shift", "difference" or "transform"
        missing="TODO: when the result is missing",
        stata=None,  # the Stata operator or command reproduced, or None for an addition
        summary="TODO: one line.",
    ),
)
def {name}(index: PanelIndex, values: FloatArray) -> FloatArray:
    """TODO: describe what is computed.

    ``values`` is ``(n_obs,)`` or ``(n_obs, k)`` in the sorted order of the panel; the result
    has the same shape and ``NaN`` where the value is missing.
    """
    # Stata: TODO command, or STATA-QUIRK: the behaviour reproduced on purpose
    raise NotImplementedError
'''

DESCRIPTION = '''"""{title} (description of a panel)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd

from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.describe._spec import DescribeSpec, description

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase


@dataclass(frozen=True, kw_only=True)
class {klass}(Result):
    """Result of :func:`{name}`."""

    table: pd.DataFrame

    def to_frame(self) -> pd.DataFrame:
        """Return the result as a table."""
        return self.table

    def summary(self, lang: Lang | None = None) -> str:
        """Return a readable report."""
        return self.table.to_string()


@description(
    "{name}",
    DescribeSpec(
        stata=None,  # the Stata command reproduced (and say "Stata" in this module), or None
        summary="TODO: one line.",
    ),
)
def {name}(panel: PanelBase) -> {klass}:
    """TODO: describe what is computed."""
    raise NotImplementedError
'''

TEST = '''"""Unit tests of {name}."""

import pandas as pd
import pytest

import pyxtpanel as xt


@pytest.fixture
def panel() -> xt.PanelData:
    frame = pd.DataFrame({{"id": [1, 1, 1, 2, 2], "t": [1, 2, 3, 1, 2], "x": [1.0, 2, 4, 3, 5]}})
    return xt.xtset(frame, "id", "t")


def test_{name}_hand_computed(panel: xt.PanelData) -> None:
    pytest.skip("TODO: a value computed by hand on this small panel")
'''

CASE = '''"""{name} against Stata."""

from __future__ import annotations

import numpy as np

from validation.registry import case

DOMAIN = "{name}"


def _python(name: str) -> dict[str, np.ndarray]:
    raise NotImplementedError


for _name in ("gaps", "unsorted_ids"):
    _tag = f"{name}_{{_name}}"

    @case(_tag, domain=DOMAIN, dataset=f"sim_{{_name}}", stata=["xtset id t"])  # TODO: Stata
    def _run(_name: str = _name) -> dict[str, np.ndarray]:
        return _python(_name)
'''

CHECK = '''"""Audit check: {title}."""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.explore.checks._context import AuditContext
from pyxtpanel.explore.checks._spec import Alert, CheckSpec, check

# one message and one advice per code, in both languages
_i18n.register(
    "explore",
    {{
        "explore.audit.{name}_problem": {{
            "en": "TODO: what was found ({{n}} units).",
            "fr": "TODO : ce qui a été trouvé ({{n}} unités).",
        }},
        "explore.audit.{name}_problem.advice": {{
            "en": "TODO: what to do about it.",
            "fr": "TODO : que faire.",
        }},
    }},
)


@check(
    "{name}",
    CheckSpec(
        summary="TODO: one line.",
        codes=("{name}_problem",),  # every code the check can raise
        order=90,  # position in the report
    ),
)
def {name}(ctx: AuditContext) -> list[Alert]:
    """TODO: describe what is checked (read the shared results from ``ctx``)."""
    return []
'''

CHECK_TEST = '''"""Unit tests of the audit check {name}."""

import pandas as pd
import pytest

import pyxtpanel as xt


@pytest.fixture
def panel() -> xt.PanelData:
    frame = pd.DataFrame({{"id": [1, 1, 1, 2, 2], "t": [1, 2, 3, 1, 2], "x": [1.0, 2, 4, 3, 5]}})
    return xt.xtset(frame, "id", "t")


def test_{name}_raises_its_alert_on_a_panel_built_for_it(panel: xt.PanelData) -> None:
    pytest.skip("TODO: build a panel that raises the alert and assert its code and level")
'''


FIGURE = '''"""``plot(kind="{name}")``: {title} (figure)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from pyxtpanel import _i18n
from pyxtpanel._types import Lang
from pyxtpanel.plot._spec import PlotSpec, plot_kind
from pyxtpanel.plot.style import chart_style, finish, label, new_figure, resolve_lang

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "plot",
    {{
        "plot.{name}.title": {{"en": "TODO: title", "fr": "TODO : titre"}},
    }},
)


@dataclass(frozen=True)
class {klass}Data:
    """What the figure draws, computed without matplotlib (so that it can be tested).

    Attributes
    ----------
    TODO
        TODO: describe the arrays the drawing reads.
    """


def {name}_data(panel: PanelBase) -> {klass}Data:
    """TODO: compute what is drawn from the panel (compare it with a table of the package)."""
    raise NotImplementedError


@plot_kind(
    "{name}",
    PlotSpec(
        summary="TODO: one line.",
        stata=None,  # the Stata command this resembles (xtline...), or None for an addition
    ),
)
def plot_{name}(panel: PanelBase, *, lang: Lang | None = None) -> Figure:
    """TODO: describe what is drawn.

    Parameters
    ----------
    panel
        The declared panel.
    lang
        ``"en"`` or ``"fr"`` (default: the global option).

    Returns
    -------
    matplotlib.figure.Figure
        The figure; save it with :func:`pyxtpanel.plot.save_figure`.
    """
    chosen = resolve_lang(lang)
    data = {name}_data(panel)
    with chart_style():
        fig, ax = new_figure()
        ax.set_title(label("plot.{name}.title", chosen))
        raise NotImplementedError(data)
        return finish(fig, chosen)
'''

FIGURE_TEST = '''"""Figure {name}: what is drawn equals what is computed."""

from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("matplotlib")

import pyxtpanel as xt  # noqa: E402


@pytest.fixture
def panel() -> xt.PanelData:
    frame = pd.DataFrame({{"id": [1, 1, 1, 2, 2], "t": [1, 2, 3, 1, 2], "x": [1.0, 2, 4, 3, 5]}})
    return xt.xtset(frame, "id", "t")


def test_{name}_data_by_hand(panel: xt.PanelData) -> None:
    pytest.skip("TODO: the computed data of this small panel, worked out by hand")


def test_{name}_draws_what_is_computed(panel: xt.PanelData) -> None:
    pytest.skip("TODO: read the artists of panel.plot({name!r}) and compare them with the data")


def test_{name}_texts_follow_the_language(panel: xt.PanelData) -> None:
    pytest.skip("TODO: the title in English and in French")
'''


DIAG_HEAD = '''"""{title} test (category {category}, family {family}).

Validated against Stata 17: TODO name the command and say which cases of
``validation/cases/diagnostics_{family}_{name}.py`` check it.
"""

from __future__ import annotations

{typing}from pyxtpanel import _i18n
from pyxtpanel.diagnostics._core import {imports}
{checking}
_i18n.register(
    "diag",
    {{
        "diag.{name}.title": {{"en": "TODO: title", "fr": "TODO : titre"}},
        "diag.{name}.null": {{"en": "TODO: null hypothesis", "fr": "TODO : hypothèse nulle"}},
        "diag.{name}.alternative": {{
            "en": "TODO: alternative hypothesis",
            "fr": "TODO : hypothèse alternative",
        }},
    }},
)


@register_test(
    "{name}",
    TestSpec(
        family="{family}",
        category="{category}",
        wave=1,
        reference="TODO: author (year), journal",
        stata="TODO: the Stata command, or 'mata' for an independent implementation",
        validation=("{name}_basic",),
        inputs="{inputs}",
        balanced="{balanced}",
        summary="TODO: one line.",
    ),
)
'''

DIAG_A = (
    DIAG_HEAD
    + '''def {name}(panel: PanelBase, variable: str) -> TestResult:
    """TODO: what the test does, its options and its validity conditions."""
    arrays = panel_variables(panel, [variable], test="{name}")
    # Stata: TODO command and options this reproduces
    raise NotImplementedError(arrays.n_obs)
'''
)

DIAG_B = (
    DIAG_HEAD
    + '''def {name}(residuals: object) -> TestResult:
    """TODO: what the test does, its arguments and its validity conditions."""
    values = as_array("residuals", residuals, ndim=1)
    # Stata: TODO command and options this reproduces
    raise NotImplementedError(values.size)
'''
)

DIAG_TEST = '''"""{title}: result, refusals and hand-computed value."""

from __future__ import annotations

import pytest


def test_{name}_hand_computed() -> None:
    pytest.skip("TODO: a value computed by hand on a tiny example")


def test_{name}_refuses_what_it_cannot_handle() -> None:
    pytest.skip("TODO: unbalanced panel, gap, too few units or non-finite arrays, as it needs")


def test_{name}_result_is_bilingual() -> None:
    pytest.skip("TODO: summary('en') and summary('fr') carry the statistic and the p-value")
'''

DIAG_CASE = '''"""{title} against Stata 17 (domain ``diagnostics_{family}``)."""

from __future__ import annotations

import numpy as np

from validation.registry import case

DOMAIN = "diagnostics_{family}"


def _python() -> dict[str, np.ndarray]:
    raise NotImplementedError("TODO: call the Python test and return the same objects")


@case("{name}_basic", domain=DOMAIN, dataset="sim_balanced", stata=["TODO: the Stata command"])
def _basic() -> dict[str, np.ndarray]:
    return _python()
'''


def create(
    kind: str,
    name: str,
    root: Path = ROOT,
    *,
    family: str | None = None,
    category: str = "A",
) -> list[Path]:
    """Write the files of a new feature (see ``KINDS``) and return their paths.

    Raises
    ------
    ValueError
        If the kind is unknown, the name is not a snake_case identifier, or a file exists.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}; available: {', '.join(KINDS)}")
    if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
        raise ValueError(f"{name!r} must be a snake_case identifier")
    if kind == "test":
        return _write(_test_files(name, family, category, root))
    package = {"operateur": "ops", "description": "describe", "controle": "explore"}.get(kind, "")
    klass = "".join(part.capitalize() for part in name.split("_"))
    title = name.replace("_", " ").capitalize()
    if kind == "graphique":
        return _write(
            {
                root / "src" / "pyxtpanel" / "plot" / "kinds" / f"{name}.py": FIGURE.format(
                    name=name, title=title, klass=klass
                ),
                root / "tests" / "plot" / f"test_figure_{name}.py": FIGURE_TEST.format(name=name),
            }
        )
    template = OPERATOR if kind == "operateur" else DESCRIPTION
    if kind == "controle":
        files = {
            root / "src" / "pyxtpanel" / "explore" / "checks" / f"{name}.py": CHECK.format(
                name=name, title=title
            ),
            root / "tests" / "unit" / f"test_explore_check_{name}.py": CHECK_TEST.format(name=name),
        }
        return _write(files)
    files = {
        root / "src" / "pyxtpanel" / package / f"{name}.py": template.format(
            name=name, title=title, klass=klass
        ),
        root / "tests" / "unit" / f"test_{package}_{name}.py": TEST.format(name=name),
        root / "validation" / "cases" / f"{name}.py": CASE.format(name=name),
    }
    return _write(files)


def _test_files(name: str, family: str | None, category: str, root: Path) -> dict[Path, str]:
    """Files of a new statistical test: module, unit test and Stata validation case."""
    from pyxtpanel.diagnostics._core.spec import FAMILIES

    if family not in FAMILIES:
        raise ValueError(f"--famille must be one of: {', '.join(FAMILIES)}")
    if category not in ("A", "B"):
        raise ValueError("--categorie must be A (on data) or B (arrays)")
    on_panel = category == "A"
    fields = {
        "name": name,
        "family": family,
        "category": category,
        "title": name.replace("_", " ").capitalize(),
        "inputs": "panel" if on_panel else "arrays",
        "balanced": "both" if on_panel else "na",
        "imports": (
            "TestResult, TestSpec, panel_variables, register_test"
            if on_panel
            else "TestResult, TestSpec, as_array, register_test"
        ),
        "typing": "from typing import TYPE_CHECKING\n\n" if on_panel else "",
        "checking": (
            "\nif TYPE_CHECKING:\n    from pyxtpanel.panel.base import PanelBase\n"
            if on_panel
            else ""
        ),
    }
    module = (DIAG_A if on_panel else DIAG_B).format(**fields)
    return {
        root / "src" / "pyxtpanel" / "diagnostics" / family / f"{name}.py": module,
        root / "tests" / "unit" / "diagnostics" / f"test_{name}.py": DIAG_TEST.format(**fields),
        root / "validation" / "cases" / f"diagnostics_{family}_{name}.py": DIAG_CASE.format(
            **fields
        ),
    }


def _write(files: dict[Path, str]) -> list[Path]:
    for path in files:
        if path.exists():
            raise ValueError(f"{path} already exists; nothing was written")
    for path, text in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")
    return list(files)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    parser.add_argument("kind", choices=list(KINDS))
    parser.add_argument("name")
    parser.add_argument("--famille", help="family of a test (folder of diagnostics/)")
    parser.add_argument("--categorie", default="A", choices=["A", "B"], help="category of a test")
    args = parser.parse_args(argv)
    try:
        written = create(args.kind, args.name, family=args.famille, category=args.categorie)
    except ValueError as error:
        print(error)
        return 1
    for path in written:
        print("created", path.relative_to(ROOT))
    if args.kind == "test":
        print("the contract test fails until the case is run in Stata and frozen:")
        print(f"  python -m validation.run diagnostics_{args.famille}")
    if args.kind == "graphique":
        print("remember: name the figure in part 6 of docs/documentation.md (table of kinds),")
        print("then python -m tools.galerie and python -m tools.doc --write")
    if args.kind == "description":
        print("remember: add a sample call to SAMPLE_CALLS in tests/contract/test_public_api.py")
    print("remember: record the feature in CHANGELOG.md, then python -m tools.stata_map --write")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
