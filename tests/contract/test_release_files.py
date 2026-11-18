"""The files that release the package agree with each other (A.1.3, A.1.5, J.1.3).

* the version is the same in ``pyproject.toml``, ``CITATION.cff`` and ``_version.py``;
* the CI matrix tests the Python versions the package announces, and the minimum versions of numpy,
  pandas and scipy that ``pyproject.toml`` promises;
* the release workflow runs the CI first and publishes through trusted publishing (no token);
* the pre-commit hooks run the same tools as the CI.

The workflow files are read as YAML, so these tests need PyYAML (a dependency of ``pre-commit``):
they are skipped where it is absent.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest

from pyxtpanel._version import __version__

yaml = pytest.importorskip("yaml")
tomllib = pytest.importorskip("tomllib")  # Python 3.11 or later: the files are checked there

ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def load(path: str) -> dict[str, Any]:
    loaded = yaml.safe_load((ROOT / path).read_text(encoding="utf-8"))
    assert isinstance(loaded, dict)
    return loaded


CI = load(".github/workflows/ci.yml")
RELEASE = load(".github/workflows/release.yml")


def test_the_version_is_the_same_everywhere() -> None:
    citation = re.search(r"^version:\s*(\S+)", (ROOT / "CITATION.cff").read_text("utf-8"), re.M)
    assert citation is not None
    assert PYPROJECT["project"]["version"] == __version__ == citation.group(1)


def test_the_ci_tests_every_python_the_package_announces() -> None:
    announced = {
        m.group(1)
        for c in PYPROJECT["project"]["classifiers"]
        if (m := re.fullmatch(r"Programming Language :: Python :: (3\.\d+)", c))
    }
    tested = set(CI["jobs"]["test"]["strategy"]["matrix"]["python"])
    assert announced == tested
    assert PYPROJECT["project"]["requires-python"] == f">={min(tested, key=_minor)}"


def _minor(version: str) -> int:
    return int(version.split(".")[1])


def test_the_ci_tests_the_minimum_versions_the_package_promises() -> None:
    minimum = {
        name: version
        for dep in PYPROJECT["project"]["dependencies"]
        if (m := re.fullmatch(r"(numpy|pandas|scipy)>=(\d+\.\d+)", dep))
        for name, version in [(m.group(1), m.group(2))]
    }
    assert set(minimum) == {"numpy", "pandas", "scipy"}
    steps = CI["jobs"]["test"]["steps"]
    install = next(s["run"] for s in steps if s.get("name", "").startswith("Install (minimum"))
    for name, version in minimum.items():
        assert f'"{name}=={version}.*"' in install, f"the minimum job must pin {name}=={version}.*"
    assert "matplotlib" not in install and "[plot]" not in install, (
        "the minimal install has no plot"
    )


def test_the_ci_has_a_job_where_linearmodels_is_installed() -> None:
    """The cross-checks with linearmodels are skipped without it: one job must install it."""
    job = CI["jobs"]["test"]
    assert any(entry["deps"] == "linearmodels" for entry in job["strategy"]["matrix"]["include"])
    install = next(s for s in job["steps"] if s.get("name") == "Install (with linearmodels)")
    assert "linearmodels" in install["run"] and install["if"] == "matrix.deps == 'linearmodels'"


def test_the_release_runs_the_ci_then_publishes_without_a_token() -> None:
    jobs = RELEASE["jobs"]
    assert jobs["checks"]["uses"] == "./.github/workflows/ci.yml"
    assert "workflow_call" in CI[True], "the CI can be called by the release workflow"
    for name in ("publish-pypi", "publish-testpypi"):
        assert jobs[name]["permissions"] == {"id-token": "write"}, "trusted publishing"
        assert jobs[name]["needs"] == "build"
    steps = yaml.safe_dump(jobs)  # the jobs only: the comments may say "no password"
    assert "secrets." not in steps and "password" not in steps, "no stored credential"
    assert jobs["publish-pypi"]["environment"]["name"] == "pypi"


def test_the_release_installs_the_published_version_from_pypi_and_tests_it() -> None:
    """After the publication a clean environment installs the tag's version and smoke-tests it."""
    job = RELEASE["jobs"]["verify-pypi"]
    assert job["needs"] == "publish-pypi" and "refs/tags" in job["if"]
    text = "\n".join(step.get("run", "") for step in job["steps"])
    assert "pyxtpanel==$version" in text, "the version of the tag, not whatever pip finds"
    assert "for attempt" in text, "the index can take minutes to list a new release"
    assert "-I smoke_install.py" in text and "GITHUB_REF_NAME#v" in text
    smoke = (ROOT / "tools" / "smoke_install.py").read_text(encoding="utf-8")
    assert "expected_version" in smoke, "the smoke test checks the installed version"


def test_the_pre_commit_hooks_run_the_tools_of_the_ci() -> None:
    config = load(".pre-commit-config.yaml")
    entries = {h["entry"] for repo in config["repos"] for h in repo["hooks"]}
    assert entries == {"ruff check", "ruff format", "mypy", "lint-imports"}
    quality = "\n".join(s.get("run", "") for s in CI["jobs"]["quality"]["steps"])
    for tool in ("ruff check", "ruff format", "mypy", "lint-imports"):
        assert tool in quality
