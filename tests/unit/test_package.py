import subprocess
import sys

import pytest

import pyxtpanel


def test_version_is_a_string() -> None:
    assert isinstance(pyxtpanel.__version__, str)


def test_import_is_light() -> None:
    code = (
        "import sys, pyxtpanel;"
        "bad = [m for m in ('matplotlib', 'numba', 'pandas', 'scipy') if m in sys.modules];"
        "print(','.join(bad))"
    )
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == ""


def test_public_names_load_on_first_use() -> None:
    for name in pyxtpanel.__all__:
        assert getattr(pyxtpanel, name) is not None
    assert "xtset" in dir(pyxtpanel)
    with pytest.raises(AttributeError, match="no attribute"):
        pyxtpanel.not_a_name  # noqa: B018


def test_all_lists_exactly_the_lazy_names_and_the_version() -> None:
    assert sorted(pyxtpanel.__all__) == sorted([*pyxtpanel._LAZY, "__version__"])
