"""Contract of the statistical tests (block L): metadata, texts and validation against Stata.

The rule that matters most is the last one: **a test without a case run in Stata 17 and frozen in
``tests/golden/`` fails here**.  There is no other way to be registered (decision of 2026-11-07).
The checks are plain functions so that the tests below can show that each one catches the
violation it exists for.
"""

from __future__ import annotations

import inspect
import json
from collections.abc import Callable, Mapping
from dataclasses import replace
from pathlib import Path

import pytest
from validation.registry import CASES, Case, discover

from pyxtpanel import _i18n
from pyxtpanel.diagnostics._core.spec import FAMILIES, TESTS, TestSpec

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden"
UNIT = ROOT / "tests" / "unit" / "diagnostics"

discover()


def golden_cases(family: str, golden_dir: Path = GOLDEN) -> tuple[str, set[str]]:
    """Return the Stata version that froze the domain of ``family`` and its frozen case tags."""
    path = golden_dir / f"diagnostics_{family}_golden.json"
    if not path.exists():
        return "", set()
    master = json.loads(path.read_text(encoding="utf-8"))
    return str(master["meta"]["stata"]), set(master["cases"])


def metadata_problems(spec: TestSpec, func: Callable[..., object]) -> list[str]:
    """Rules on the declaration: family folder, category, inputs, non-empty fields."""
    found: list[str] = []
    module = inspect.getmodule(func)
    if spec.family not in FAMILIES:
        found.append(f"unknown family {spec.family!r}")
    elif module is None or not module.__name__.startswith(f"pyxtpanel.diagnostics.{spec.family}."):
        found.append(f"must live in diagnostics/{spec.family}/")
    if (spec.category, spec.inputs) not in (("A", "panel"), ("B", "arrays")):
        found.append("category A takes a panel, category B takes arrays")
    if (spec.inputs == "arrays") != (spec.balanced == "na"):
        found.append("balanced is 'na' for arrays and 'required' or 'both' for a panel")
    found += [
        f"{f} is empty" for f in ("reference", "stata", "summary") if not getattr(spec, f).strip()
    ]
    return found


def text_problems(name: str, func: Callable[..., object], unit_dir: Path) -> list[str]:
    """Rules on documentation: docstring, messages in both languages, provenance, unit test."""
    found: list[str] = []
    if not inspect.getdoc(func):
        found.append("no docstring")
    for suffix in ("title", "null", "alternative"):
        for lang in ("en", "fr"):
            key = f"diag.{name}.{suffix}"
            if key not in _i18n.all_keys() or not _i18n.has_language(key, lang):  # type: ignore[arg-type]
                found.append(f"message diag.{name}.{suffix} missing in {lang}")
    module = inspect.getmodule(func)
    source = inspect.getsource(module) if module is not None else ""
    if "# Stata:" not in source:
        found.append("the code carries no '# Stata:' provenance comment")
    if "validated" not in source.lower():
        found.append("the module does not say how it is validated")
    if not (unit_dir / f"test_{name}.py").exists():
        found.append(f"tests/unit/diagnostics/test_{name}.py is missing")
    return found


def validation_problems(spec: TestSpec, cases: Mapping[str, Case], golden_dir: Path) -> list[str]:
    """Validation against Stata: not negotiable."""
    found: list[str] = []
    if not spec.validation:
        found.append("no validation case: every test is validated against Stata 17")
    version, frozen = golden_cases(spec.family, golden_dir)
    if spec.validation and not version.startswith("Stata 17"):
        found.append(f"no golden file diagnostics_{spec.family}_golden.json produced by Stata 17")
    commands = [part.strip() for part in spec.stata.split(";") if part.strip()]
    run_by_cases: set[str] = set()
    for tag in spec.validation:
        case = cases.get(tag)
        if case is None:
            found.append(f"validation case {tag!r} does not exist")
            continue
        if case.domain != f"diagnostics_{spec.family}":
            found.append(f"case {tag!r} must be in the domain diagnostics_{spec.family}")
        text = " ".join(case.stata)
        used = {c for c in commands if c in text}
        run_by_cases |= used
        if commands and not used:
            found.append(f"case {tag!r} never runs the Stata command {commands[0]!r}")
        if tag not in frozen:
            found.append(f"case {tag!r} is not frozen from Stata (run python -m validation.run)")
    if spec.validation:
        found += [
            f"no case runs the Stata command {c!r}" for c in commands if c not in run_by_cases
        ]
    return found


def problems(
    name: str,
    spec: TestSpec,
    func: Callable[..., object],
    cases: Mapping[str, Case] = CASES,
    golden_dir: Path = GOLDEN,
    unit_dir: Path = UNIT,
) -> list[str]:
    """Return every rule the registered test ``name`` breaks (empty when it is complete)."""
    return [
        *metadata_problems(spec, func),
        *text_problems(name, func, unit_dir),
        *validation_problems(spec, cases, golden_dir),
    ]


@pytest.mark.parametrize("name", TESTS.names())
def test_each_registered_test_is_complete_and_validated_in_stata(name: str) -> None:
    assert problems(name, TESTS.spec(name), TESTS.get(name)) == []


def test_the_registry_is_empty_or_every_test_is_in_a_known_family() -> None:
    assert all(spec.family in FAMILIES for _, spec in TESTS.specs())


# ---------------------------------------------------------------- the rules catch their violation

SPEC = TestSpec(
    family="unitroot",
    category="A",
    wave=1,
    reference="Author (2000)",
    stata="xtunitroot ips",
    validation=("ips_balanced",),
    inputs="panel",
    summary="a test",
    balanced="both",
)


def probe() -> None:
    """A probe test."""


def case_for(*lines: str, domain: str = "diagnostics_unitroot") -> Case:
    return Case("ips_balanced", domain, "x", tuple(lines), lambda: {})


@pytest.fixture
def golden_dir(tmp_path: Path) -> Path:
    (tmp_path / "diagnostics_unitroot_golden.json").write_text(
        json.dumps({"meta": {"stata": "Stata 17.0 MP"}, "cases": {"ips_balanced": {}}}),
        encoding="utf-8",
    )
    return tmp_path


def run(spec: TestSpec, cases: Mapping[str, Case], golden_dir: Path) -> list[str]:
    return problems("probe", spec, probe, cases, golden_dir, golden_dir)


def test_a_spec_without_validation_case_is_refused(golden_dir: Path) -> None:
    found = run(replace(SPEC, validation=()), {}, golden_dir)
    assert any("every test is validated against Stata 17" in line for line in found)


def test_an_unknown_validation_case_is_refused(golden_dir: Path) -> None:
    assert any("does not exist" in line for line in run(SPEC, {}, golden_dir))


def test_a_declared_stata_command_that_no_case_runs_is_refused(golden_dir: Path) -> None:
    two = replace(SPEC, stata="xtunitroot ips; xtabond2")
    found = run(two, {"ips_balanced": case_for("xtunitroot ips y")}, golden_dir)
    assert any("no case runs the Stata command 'xtabond2'" in line for line in found)
    ok = run(
        replace(SPEC, stata="xtunitroot ips; xtabond2"),
        {"ips_balanced": case_for("xtunitroot ips y", "xtabond2 y")},
        golden_dir,
    )
    assert not any("no case runs" in line for line in ok)


def test_a_case_that_never_runs_the_stata_command_is_refused(golden_dir: Path) -> None:
    found = run(SPEC, {"ips_balanced": case_for("regress y x")}, golden_dir)
    assert any("never runs the Stata command 'xtunitroot ips'" in line for line in found)


def test_a_case_in_the_wrong_domain_is_refused(golden_dir: Path) -> None:
    cases = {"ips_balanced": case_for("xtunitroot ips y", domain="explore")}
    assert any("must be in the domain" in line for line in run(SPEC, cases, golden_dir))


def test_a_case_not_frozen_from_stata_is_refused(golden_dir: Path) -> None:
    other = replace(SPEC, validation=("ips_balanced", "ips_gaps"))
    cases = {
        "ips_balanced": case_for("xtunitroot ips y"),
        "ips_gaps": Case("ips_gaps", "diagnostics_unitroot", "x", ("xtunitroot ips y",), dict),
    }
    found = run(other, cases, golden_dir)
    assert any("'ips_gaps' is not frozen" in line for line in found)


def test_a_golden_file_not_produced_by_stata_17_is_refused(tmp_path: Path) -> None:
    (tmp_path / "diagnostics_unitroot_golden.json").write_text(
        json.dumps({"meta": {"stata": "python"}, "cases": {"ips_balanced": {}}}), encoding="utf-8"
    )
    cases = {"ips_balanced": case_for("xtunitroot ips y")}
    assert any("produced by Stata 17" in line for line in run(SPEC, cases, tmp_path))


def test_a_missing_golden_file_is_refused(tmp_path: Path) -> None:
    cases = {"ips_balanced": case_for("xtunitroot ips y")}
    assert any("produced by Stata 17" in line for line in run(SPEC, cases, tmp_path))


def test_an_independent_mata_reference_is_accepted_when_the_case_runs_mata(
    golden_dir: Path,
) -> None:
    mata = replace(SPEC, stata="mata")
    cases = {"ips_balanced": case_for("mata: pyx_export(...)")}
    found = run(mata, cases, golden_dir)
    assert not any("validation" in line.lower() or "Stata" in line for line in found)


def test_incoherent_metadata_is_refused(golden_dir: Path) -> None:
    assert any(
        "category A takes a panel" in line
        for line in run(replace(SPEC, category="B"), {}, golden_dir)
    )
    assert any(
        "balanced is 'na'" in line for line in run(replace(SPEC, balanced="na"), {}, golden_dir)
    )
