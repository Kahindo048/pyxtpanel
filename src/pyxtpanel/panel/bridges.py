"""Bridges to the ``linearmodels`` panel container, by duck typing (no hard dependency).

``linearmodels.panel.data.PanelData`` is a different class from :class:`pyxtpanel.PanelData`
that happens to share its name.  Both wrap a DataFrame indexed by (entity, time).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pandas as pd

from pyxtpanel._errors import SpecificationError

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase


def frame_of_linearmodels(obj: Any) -> pd.DataFrame:
    """Return the (entity, time)-indexed DataFrame held by a ``linearmodels`` panel container.

    Anything with a two-level ``dataframe`` attribute is accepted, as is a DataFrame itself.

    Raises
    ------
    SpecificationError
        If ``obj`` holds no DataFrame with a two-level index.
    """
    frame = obj if isinstance(obj, pd.DataFrame) else getattr(obj, "dataframe", None)
    if not isinstance(frame, pd.DataFrame) or frame.index.nlevels != 2:
        raise SpecificationError(
            "expected a linearmodels PanelData (or a DataFrame) indexed by (entity, time)",
            stata_code=198,
            key="bridges.not_linearmodels",
        )
    return frame


def to_linearmodels(panel: PanelBase) -> Any:
    """Return the data as a ``linearmodels.panel.data.PanelData``.

    Raises
    ------
    ImportError
        If ``linearmodels`` is not installed.
    """
    try:
        from linearmodels.panel.data import PanelData as LinearModelsPanelData
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise ImportError("linearmodels is not installed: pip install linearmodels") from error
    return LinearModelsPanelData(panel.to_frame(order="sorted", index="keys"))
