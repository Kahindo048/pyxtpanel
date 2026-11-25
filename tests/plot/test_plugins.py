"""Figure kinds supplied by other packages through the entry-point group (B.1.3).

No real package is installed: the discovery is replaced by stand-ins that have what an
``importlib.metadata.EntryPoint`` has (``name``, ``value``, ``load()``).
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from typing import Any

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")

import pyxtpanel as xt  # noqa: E402
from pyxtpanel._errors import PyxtpanelError, SpecificationError  # noqa: E402
from pyxtpanel.plot import KINDS, _plugins  # noqa: E402


class Source:
    """What an entry point of the group points to: ``names()`` and ``get(name)``."""

    def __init__(self, kinds: dict[str, Callable[..., Any]]) -> None:
        self.kinds = kinds
        self.loaded: list[str] = []

    def names(self) -> list[str]:
        return list(self.kinds)

    def get(self, name: str) -> Callable[..., Any]:
        self.loaded.append(name)
        return self.kinds[name]


class Entry:
    def __init__(self, name: str, target: object | Exception) -> None:
        self.name, self.value, self._target = name, f"fake.{name}:KINDS", target
        self.loads = 0

    def load(self) -> object:
        self.loads += 1
        if isinstance(self._target, Exception):
            raise self._target
        return self._target


def install(monkeypatch: pytest.MonkeyPatch, *entries: Entry) -> None:
    monkeypatch.setattr(_plugins, "_discover", lambda: sorted(entries, key=lambda e: e.name))
    _plugins._reset()


@pytest.fixture(autouse=True)
def fresh_discovery() -> Iterator[None]:
    _plugins._reset()
    yield
    _plugins._reset()


@pytest.fixture(scope="module")
def panel() -> xt.PanelData:
    frame = pd.DataFrame(
        {"id": np.repeat([1, 2], 3), "t": np.tile([1, 2, 3], 2), "x": np.arange(6.0)}
    )
    return xt.xtset(frame, "id", "t")


def marker(panel: object, variable: str, *, bins: int = 10, lang: str | None = None) -> str:
    return f"drawn {variable} {bins}"


def test_a_plugin_kind_is_drawn_through_plot_and_the_method(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    install(monkeypatch, Entry("viz", Source({"marker": marker})))
    assert xt.plot.plot(panel, "marker", variable="x", bins=3) == "drawn x 3"
    assert panel.plot("MARKER", variable="x") == "drawn x 10"  # case does not matter


def test_a_plugin_call_is_checked_against_its_signature(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    install(monkeypatch, Entry("viz", Source({"marker": marker})))
    with pytest.raises(SpecificationError) as error:
        panel.plot("marker", variable="x", bogus=1)
    assert error.value.key == "plot.bad_call"
    assert "bogus" in str(error.value) and "bins=10" in str(error.value)


def test_a_builtin_kind_never_touches_the_plugins(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    entry = Entry("viz", Source({"lines": marker}))
    install(monkeypatch, entry)
    figure = panel.plot("lines", variable="x")
    assert entry.loads == 0  # lazy: nothing was looked up
    assert figure.__class__.__name__ == "Figure"


def test_the_builtin_kind_wins_over_a_plugin_of_the_same_name(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    source = Source({"lines": marker})
    install(monkeypatch, Entry("viz", source))
    assert panel.plot("lines", variable="x").__class__.__name__ == "Figure"
    assert source.loaded == []


def test_among_plugins_the_first_entry_name_wins(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    def other(panel: object, variable: str) -> str:
        return "other"

    install(
        monkeypatch,
        Entry("b_second", Source({"marker": other})),
        Entry("a_first", Source({"marker": marker})),
    )
    assert panel.plot("marker", variable="x") == "drawn x 10"


def test_a_broken_plugin_does_not_break_the_builtin_figures(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    install(monkeypatch, Entry("broken", ImportError("no module named foo")))
    assert panel.plot("lines", variable="x").__class__.__name__ == "Figure"


@pytest.mark.parametrize("lang", ["en", "fr"])
def test_a_kind_that_may_be_in_a_broken_plugin_gives_a_clear_error(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData, lang: str
) -> None:
    install(monkeypatch, Entry("broken", ImportError("no module named foo")))
    with pytest.raises(PyxtpanelError) as error:
        panel.plot("marker", variable="x", lang=lang)
    assert not isinstance(error.value, ImportError)
    assert error.value.key == "plot.plugin_broken"
    text = str(error.value)
    assert "broken" in text and "fake.broken:KINDS" in text and "no module named foo" in text
    assert ("inutilisable" if lang == "fr" else "cannot be used") in text


def test_a_working_plugin_is_used_even_if_another_one_is_broken(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    install(
        monkeypatch,
        Entry("a_broken", RuntimeError("boom")),
        Entry("b_good", Source({"marker": marker})),
    )
    assert panel.plot("marker", variable="x") == "drawn x 10"


def test_a_plugin_whose_get_fails_is_reported(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    class Failing(Source):
        def get(self, name: str) -> Callable[..., Any]:
            raise KeyError(name)

    install(monkeypatch, Entry("viz", Failing({"marker": marker})))
    with pytest.raises(PyxtpanelError) as error:
        panel.plot("marker", variable="x")
    assert error.value.key == "plot.plugin_broken"


def test_an_object_without_names_is_a_broken_plugin(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    install(monkeypatch, Entry("viz", object()))
    with pytest.raises(PyxtpanelError) as error:
        panel.plot("marker", variable="x")
    assert error.value.key == "plot.plugin_broken"


def test_an_unknown_kind_names_pyxtviz_when_no_package_supplies_kinds(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    install(monkeypatch)
    with pytest.raises(SpecificationError) as error:
        panel.plot("nothing", variable="x")
    assert error.value.key == "plot.unknown_kind_install"
    text = str(error.value)
    assert "pyxtviz" in text and "pip install pyxtviz" in text
    assert all(name in text for name in KINDS.names())


def test_an_unknown_kind_lists_the_plugin_kinds_too(
    monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData
) -> None:
    install(monkeypatch, Entry("viz", Source({"marker": marker})))
    with pytest.raises(SpecificationError) as error:
        panel.plot("nothing", variable="x")
    assert error.value.key == "plot.unknown_kind"
    assert "marker" in str(error.value) and "lines" in str(error.value)


def test_the_discovery_is_done_once(monkeypatch: pytest.MonkeyPatch, panel: xt.PanelData) -> None:
    calls = []
    entry = Entry("viz", Source({"marker": marker}))
    monkeypatch.setattr(_plugins, "_discover", lambda: calls.append(1) or [entry])
    _plugins._reset()
    panel.plot("marker", variable="x")
    panel.plot("marker", variable="x")
    assert len(calls) == 1


def test_the_real_discovery_runs_without_a_plugin_installed() -> None:
    assert isinstance(_plugins._discover(), list)
