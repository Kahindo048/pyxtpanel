"""Invariance properties of the declaration (C.6.4)."""

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxtpanel._errors import PyxtpanelWarning
from pyxtpanel.panel.base import PanelBase

pairs = st.lists(
    st.tuples(st.integers(0, 6), st.integers(0, 12)), min_size=1, max_size=45, unique=True
)
quiet = pytest.mark.filterwarnings("ignore::pyxtpanel._errors.PyxtpanelWarning")


def declare(data, **kw):  # type: ignore[no-untyped-def]
    frame = pd.DataFrame(data, columns=["id", "t"])
    return PanelBase.from_frame(frame, "id", "t", **kw)


@quiet
@settings(max_examples=60, deadline=None)
@given(data=pairs, shift=st.integers(-1000, 1000))
def test_shifting_the_time_does_not_change_the_structure(data, shift) -> None:  # type: ignore[no-untyped-def]
    a, b = declare(data), declare([(i, t + shift) for i, t in data])
    assert (a.index.rank == b.index.rank).all() and (a.index.offsets == b.index.offsets).all()
    assert (a.balanced, a.n_gaps) == (b.balanced, b.n_gaps)
    assert b.t_min == a.t_min + shift


@quiet
@settings(max_examples=60, deadline=None)
@given(data=pairs)
def test_text_and_integer_identifiers_give_the_same_structure(data) -> None:  # type: ignore[no-untyped-def]
    ints = declare(data)
    text = declare([(f"unit{i:03d}", t) for i, t in data])
    assert (ints.index.offsets == text.index.offsets).all()
    assert (ints.index.rank == text.index.rank).all()
    assert text.unit_labels.tolist() == [f"unit{i:03d}" for i in ints.unit_labels]


@quiet
@settings(max_examples=60, deadline=None)
@given(data=pairs, factor=st.integers(2, 5))
def test_scaling_the_time_with_delta_does_not_change_the_structure(data, factor) -> None:  # type: ignore[no-untyped-def]
    a = declare(data)
    b = declare([(i, t * factor) for i, t in data], delta=factor)
    assert (a.index.rank == b.index.rank).all() and a.n_gaps == b.n_gaps


@quiet
@settings(max_examples=40, deadline=None)
@given(data=pairs)
def test_tsfill_full_always_gives_a_strongly_balanced_panel(data) -> None:  # type: ignore[no-untyped-def]
    filled = declare(data).tsfill(full=True)
    assert filled.balanced == "strongly balanced" and not filled.has_gaps
    assert filled.n_obs == filled.n_units * filled.n_periods


@quiet
@settings(max_examples=40, deadline=None)
@given(data=pairs)
def test_selecting_everything_changes_nothing(data) -> None:  # type: ignore[no-untyped-def]
    panel = declare(data)
    same = panel.select(units=list(panel.unit_labels))
    assert (same.index.rank == panel.index.rank).all() and same.n_obs == panel.n_obs
    assert np.array_equal(same.index.source_row, panel.index.source_row)


@given(
    garbage=st.sampled_from(
        [
            pd.DataFrame({"id": [1, 1], "t": [np.nan, np.nan]}),
            pd.DataFrame({"id": [np.nan, 1], "t": [1, 1]}),
            pd.DataFrame({"id": [1, 1], "t": [1, 1]}),
            pd.DataFrame({"id": [1, 1], "t": [1.5, 2]}),
            pd.DataFrame({"id": [1, 1], "t": ["a", "b"]}),
            pd.DataFrame({"id": [], "t": []}),
        ]
    )
)
def test_damaged_input_is_refused_with_a_package_error(garbage) -> None:  # type: ignore[no-untyped-def]
    from pyxtpanel._errors import PyxtpanelError

    with pytest.raises(PyxtpanelError):
        PanelBase.from_frame(garbage, "id", "t")


def test_warning_class_is_exported() -> None:
    assert issubclass(PyxtpanelWarning, UserWarning)
