"""The generator of new features writes valid, registrable modules without touching others."""

import ast

import pytest
from tools import nouveau


@pytest.mark.parametrize("kind", ["operateur", "description"])
def test_generated_files_are_valid_python_and_registered(kind: str, tmp_path) -> None:  # type: ignore[no-untyped-def]
    written = nouveau.create(kind, "moving_thing", tmp_path)
    assert len(written) == 3 and all(p.exists() for p in written)
    for path in written:
        ast.parse(path.read_text(encoding="utf-8"))
    module = next(p for p in written if p.parent.name in ("ops", "describe"))
    text = module.read_text(encoding="utf-8")
    assert ("@operator(" if kind == "operateur" else "@description(") in text
    assert '"moving_thing"' in text


def test_the_generator_refuses_bad_names_unknown_kinds_and_existing_files(tmp_path) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(ValueError, match="snake_case"):
        nouveau.create("operateur", "Bad-Name", tmp_path)
    with pytest.raises(ValueError, match="unknown kind"):
        nouveau.create("bogus", "x", tmp_path)
    nouveau.create("operateur", "again", tmp_path)
    with pytest.raises(ValueError, match="already exists"):
        nouveau.create("operateur", "again", tmp_path)


def test_a_new_audit_check_gets_a_module_and_a_test(tmp_path) -> None:  # type: ignore[no-untyped-def]
    written = nouveau.create("controle", "short_units", tmp_path)
    assert [p.parent.name for p in written] == ["checks", "unit"]
    text = written[0].read_text(encoding="utf-8")
    assert (
        "@check(" in text
        and '"short_units"' in text
        and "explore.audit.short_units_problem" in text
    )
    for path in written:
        ast.parse(path.read_text(encoding="utf-8"))


@pytest.mark.parametrize("category", ["A", "B"])
def test_a_new_statistical_test_gets_module_unit_test_and_stata_case(  # type: ignore[no-untyped-def]
    category: str, tmp_path
) -> None:
    written = nouveau.create("test", "demo", tmp_path, family="effects", category=category)
    assert [p.name for p in written] == [
        "demo.py",
        "test_demo.py",
        "diagnostics_effects_demo.py",
    ]
    module, _, case = (p.read_text(encoding="utf-8") for p in written)
    for path in written:
        ast.parse(path.read_text(encoding="utf-8"))
    assert "@register_test(" in module and 'family="effects"' in module
    assert f'category="{category}"' in module and "# Stata:" in module
    assert ("panel_variables" if category == "A" else "as_array") in module
    assert 'DOMAIN = "diagnostics_effects"' in case and '@case("demo_basic"' in case


def test_a_new_statistical_test_needs_a_known_family_and_category(tmp_path) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(ValueError, match="--famille must be one of"):
        nouveau.create("test", "demo", tmp_path, family="nowhere")
    with pytest.raises(ValueError, match="--categorie must be A"):
        nouveau.create("test", "demo", tmp_path, family="effects", category="C")
    assert not list(tmp_path.rglob("*.py"))
