"""Public entry point for the packages that extend ``pyxtpanel`` (plots, estimators, reports).

A companion package (``pyxtviz``, an estimator...) builds on the same pieces as the foundation
itself.  This module is where they are published, so that a companion never imports a private
module (one whose name starts with an underscore): those may change in any release.

* :class:`Registry`: named entries with metadata, for a package's own extension points.
* :func:`register_messages`, :func:`message`: bilingual (English and French) messages, one table
  per domain prefix.  Use your own prefix (``viz.``, ``abond.``).
* :func:`call_reason`, :func:`call_arguments`: the wording of a call that does not fit the
  signature of the function it reaches.
* :data:`Lang`, :data:`LANGUAGES`, :data:`FloatArray`, :data:`IntArray`, :data:`BoolArray`: types.
* :class:`ExtendedPanelLike`: what a package may rely on in a panel, calendar and terms included.

The panel itself is read through :meth:`pyxtpanel.PanelData.arrays` and
:meth:`pyxtpanel.PanelData.term_arrays` (see :mod:`pyxtpanel.contract`), and its calendar through
:attr:`pyxtpanel.PanelData.calendar`.

:data:`EXTEND_VERSION` is the version of this surface.  Adding a name does not change it; removing
or changing one does.  A companion package checks it when it is imported.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

from pyxtpanel import _i18n, _signature
from pyxtpanel._registry import Registry
from pyxtpanel._types import LANGUAGES, BoolArray, FloatArray, IntArray, Lang

if TYPE_CHECKING:
    import numpy as np

    from pyxtpanel.panel.arrays import PanelArrays
    from pyxtpanel.time.calendar import Calendar

__all__ = [
    "EXTEND_VERSION",
    "LANGUAGES",
    "BoolArray",
    "ExtendedPanelLike",
    "FloatArray",
    "IntArray",
    "Lang",
    "Registry",
    "call_arguments",
    "call_reason",
    "has_language",
    "message",
    "message_fields",
    "message_keys",
    "register_messages",
]

#: Version of the extension surface; a companion package checks it when it is imported.
EXTEND_VERSION = 1

# The same objects as the internal ones (a test checks it), published under stable names.
register_messages = _i18n.register
message = _i18n.message
message_keys = _i18n.all_keys
has_language = _i18n.has_language
message_fields = _i18n.fields_of
call_reason = _signature.reason
call_arguments = _signature.accepted


@runtime_checkable
class ExtendedPanelLike(Protocol):
    """What a companion package may rely on in a declared panel.

    Every :class:`pyxtpanel.PanelData` satisfies it.  It is :class:`pyxtpanel.PanelLike` (names,
    sizes, :meth:`arrays`) plus the calendar, the list of numeric variables, the reader of terms
    and the reader of columns of any type.
    """

    id_name: str
    time_name: str

    @property
    def n_units(self) -> int:
        """Number of units."""
        ...

    @property
    def n_obs(self) -> int:
        """Number of observations."""
        ...

    @property
    def calendar(self) -> Calendar:
        """Read-only view of the calendar shared by all units."""
        ...

    @property
    def numeric_names(self) -> tuple[str, ...]:
        """Numeric variables other than the two keys, in column order (what ``arrays`` can read)."""
        ...

    def arrays(self, columns: Any, *, where: Any = None) -> PanelArrays:
        """Return the array contract for variables (see :class:`pyxtpanel.contract.PanelArrays`)."""
        ...

    def sorted_column(self, name: str) -> np.ndarray:
        """Return a column of any type (a label, a date) in the sorted order of the rows."""
        ...

    def term_arrays(self, items: Any) -> PanelArrays:
        """Return the same contract for variables or Stata-style terms (``"L.x"``, ``"D.x"``)."""
        ...
