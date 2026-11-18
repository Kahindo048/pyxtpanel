"""Random scenarios against Stata (J.4.1): one battery of commands on 160 hostile panels.

Stata runs ``pyx_battery`` on each scenario of ``sim_campaign.dta`` (declaration, ``xtsum``,
six operators, within and between, ``xtdescribe``, ``tsfill``) and exports one row each; the
package computes the same row; every cell must agree.  See ``validation/campaign.py``.
"""

from __future__ import annotations

import numpy as np

from validation import campaign
from validation.registry import case, export, register_prelude

DOMAIN = "campaign"

register_prelude(DOMAIN, campaign.MATA + "\n" + campaign.stata_program())

_LINES = [
    f"mata: R = J(0, {len(campaign.COLUMNS)}, .)",
    *[f'pyx_battery {s.k} {s.delta} "{s.fmt}"' for s in campaign.load()],
    export("campaign", "battery", "R"),
]


@case("campaign", domain=DOMAIN, dataset="sim_campaign", stata=_LINES, atol=1e-9, rtol=1e-9)
def _battery() -> dict[str, np.ndarray]:
    rows = [campaign.python_row(s.frame, s.delta, s.fmt) for s in campaign.load()]
    return {"battery": np.vstack(rows)}
