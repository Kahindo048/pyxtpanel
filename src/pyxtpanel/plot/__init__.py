"""Optional plotting extra (matplotlib is imported only here, and only when a figure is drawn).

Layer 5.  ``plot(panel, kind, ...)`` draws one of the registered kinds (``KINDS.names()``);
``PanelData.plot`` is the same call.  The look is the one of ``pyxtdpthresh``
(:mod:`pyxtpanel.plot.style`).  Importing this package does not import matplotlib.
"""

from __future__ import annotations

from pyxtpanel.plot._spec import KINDS, PlotSpec, plot, plot_kind
from pyxtpanel.plot.style import save_figure

__all__ = ["KINDS", "PlotSpec", "plot", "plot_kind", "save_figure"]
