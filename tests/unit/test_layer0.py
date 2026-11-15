import sys
from dataclasses import FrozenInstanceError, dataclass

import pandas as pd
import pytest

from pyxtpanel import _i18n
from pyxtpanel._config import get_options, option_context, set_options
from pyxtpanel._errors import (
    DataError,
    DeclarationError,
    NumericalError,
    PyxtpanelError,
    PyxtpanelWarning,
    SpecificationError,
    TimeAxisError,
)
from pyxtpanel._registry import Registry
from pyxtpanel._results import Result

# ---------------------------------------------------------------- errors


def test_error_hierarchy_and_builtin_bases() -> None:
    assert issubclass(DeclarationError, DataError) and issubclass(TimeAxisError, DataError)
    assert issubclass(DataError, PyxtpanelError) and issubclass(DataError, ValueError)
    assert issubclass(SpecificationError, ValueError)
    assert issubclass(NumericalError, ArithmeticError)
    assert issubclass(PyxtpanelWarning, UserWarning)


def test_error_carries_stata_code_and_key() -> None:
    err = DeclarationError("repeated pair", stata_code=451, key="panel.duplicates")
    assert (str(err), err.stata_code, err.key) == ("repeated pair", 451, "panel.duplicates")
    assert SpecificationError("bad option").stata_code == 198
    assert DataError("x").stata_code is None


def test_a_foundation_error_is_caught_by_an_estimator_root() -> None:
    class EstimatorError(PyxtpanelError):
        pass

    with pytest.raises(PyxtpanelError):
        raise DataError("from the foundation")
    assert not issubclass(EstimatorError, DataError)


# ---------------------------------------------------------------- options


def test_options_default_and_restore() -> None:
    assert get_options().lang == "en"
    with option_context(lang="fr") as options:
        assert options.lang == "fr"
    assert get_options().lang == "en"


def test_options_restored_after_an_exception() -> None:
    with pytest.raises(RuntimeError), option_context(lang="fr"):
        raise RuntimeError
    assert get_options().lang == "en"


def test_set_options_rejects_unknown_name_and_value() -> None:
    with pytest.raises(ValueError, match="Unknown option"):
        set_options(colour="blue")
    with pytest.raises(ValueError, match="lang must be"):
        set_options(lang="de")
    assert get_options().lang == "en"


# ---------------------------------------------------------------- messages


def test_message_languages_and_fallback() -> None:
    _i18n.register("test_l0", {"test_l0.hello": {"en": "Hello {who}", "fr": "Bonjour {who}"}})
    _i18n.register("test_l0", {"test_l0.only_en": {"en": "Only English"}})
    assert _i18n.message("test_l0.hello", "fr", who="Ada") == "Bonjour Ada"
    assert _i18n.message("test_l0.hello", who="Ada") == "Hello Ada"
    with option_context(lang="fr"):
        assert _i18n.message("test_l0.hello", who="Ada") == "Bonjour Ada"
        assert _i18n.message("test_l0.only_en") == "Only English"
    assert _i18n.has_language("test_l0.hello", "fr")
    assert not _i18n.has_language("test_l0.only_en", "fr")
    assert _i18n.fields_of("test_l0.hello", "fr") == {"who"}


def test_register_rejects_wrong_prefix_and_conflicting_duplicate() -> None:
    with pytest.raises(ValueError, match="must start with"):
        _i18n.register("test_l0", {"other.key": {"en": "x", "fr": "x"}})
    _i18n.register("test_l0", {"test_l0.dup": {"en": "a", "fr": "a"}})
    _i18n.register("test_l0", {"test_l0.dup": {"en": "a", "fr": "a"}})  # same text: harmless
    with pytest.raises(ValueError, match="already registered"):
        _i18n.register("test_l0", {"test_l0.dup": {"en": "b", "fr": "b"}})


def test_unknown_message_key_raises() -> None:
    with pytest.raises(KeyError):
        _i18n.message("test_l0.absent")


# ---------------------------------------------------------------- registry


@dataclass(frozen=True)
class _Spec:
    family: str


def test_registry_register_get_names_spec() -> None:
    reg: Registry[_Spec] = Registry("thing")

    @reg.register("Alpha", _Spec("a"))
    def alpha() -> int:
        return 1

    @reg.register("beta", _Spec("b"))
    def beta() -> int:
        return 2

    assert reg.names() == ["alpha", "beta"] and len(reg) == 2
    assert "ALPHA" in reg and "gamma" not in reg
    assert reg.get("alpha")() == 1 and reg.spec("Beta") == _Spec("b")
    assert [n for n, _ in reg.specs()] == ["alpha", "beta"]


def test_registry_unknown_and_duplicate() -> None:
    reg: Registry[_Spec] = Registry("thing")

    @reg.register("alpha", _Spec("a"))
    def alpha() -> None: ...

    with pytest.raises(SpecificationError, match="available: alpha") as info:
        reg.get("zeta")
    assert info.value.key == "registry.unknown"

    with pytest.raises(SpecificationError, match="already registered"):

        @reg.register("alpha", _Spec("a"))
        def another() -> None: ...


def test_registry_load_all_imports_submodules(tmp_path, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    pkg = tmp_path / "fakefeatures"
    (pkg / "sub").mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "sub" / "__init__.py").write_text("")
    (pkg / "one.py").write_text(
        "from fakeregistry import REG, Spec\n@REG.register('one', Spec('x'))\ndef one(): return 1\n"
    )
    (pkg / "sub" / "two.py").write_text(
        "from fakeregistry import REG, Spec\n@REG.register('two', Spec('y'))\ndef two(): return 2\n"
    )
    (tmp_path / "fakeregistry.py").write_text(
        "from dataclasses import dataclass\n"
        "from pyxtpanel._registry import Registry\n"
        "@dataclass(frozen=True)\n"
        "class Spec:\n    family: str\n"
        "REG = Registry('thing', 'fakefeatures')\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    import fakeregistry

    assert fakeregistry.REG.names() == ["one", "two"]
    for name in ("fakefeatures", "fakefeatures.one", "fakefeatures.sub.two", "fakeregistry"):
        sys.modules.pop(name, None)


# ---------------------------------------------------------------- result base


@dataclass(frozen=True, kw_only=True)
class _R(Result):
    value: float = 1.5

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame({"value": [self.value]})

    def summary(self, lang=None) -> str:  # type: ignore[no-untyped-def]
        return f"value = {self.value}"


def test_result_is_frozen_and_renders() -> None:
    r = _R()
    assert str(r) == "value = 1.5" and "1.5" in r._repr_html_()
    with pytest.raises(FrozenInstanceError):
        r.value = 2.0  # type: ignore[misc]
    with pytest.raises(NotImplementedError):
        Result().to_frame()
    with pytest.raises(TypeError):
        _R(1.5)  # type: ignore[misc]  # fields are keyword-only
