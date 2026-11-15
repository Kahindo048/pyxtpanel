"""Public entry point of the array contract used by the estimation packages.

An estimator checks ``CONTRACT_VERSION`` when it is imported and receives ``PanelArrays`` from
:meth:`pyxtpanel.PanelData.arrays`.  The guarantees are documented in
:mod:`pyxtpanel.panel.arrays`.
"""

from __future__ import annotations

from pyxtpanel.panel.arrays import CONTRACT_VERSION, PanelArrays

__all__ = ["CONTRACT_VERSION", "PanelArrays"]
