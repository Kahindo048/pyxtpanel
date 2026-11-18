"""``xttrans``: transition matrix of a variable between consecutive observations (E.3.2).

Stata pairs each observation with the **next row of the same unit**, so a gap in the calendar
is bridged (a unit observed in 1, 2 and 4 contributes the transitions 1→2 and 2→4).  A row with
a missing value is not bridged: the transitions that touch it are dropped.  ``gaps="strict"``
keeps only transitions between consecutive calendar periods; the number of bridged transitions
is always reported.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import FloatArray, Lang
from pyxtpanel.describe._format import number
from pyxtpanel.describe._spec import DescribeSpec, description
from pyxtpanel.ops.frame import term_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "describe",
    {
        "describe.xr.no_transition": {
            "en": "{name!r} has no transition between two observed values.",
            "fr": "« {name} » n'a aucune transition entre deux valeurs observées.",
        },
        "describe.xr.bad_gaps": {
            "en": "gaps must be 'bridge' or 'strict', got {gaps!r}.",
            "fr": "gaps doit valoir 'bridge' ou 'strict', reçu {gaps!r}.",
        },
        "describe.xr.bad_bins": {
            "en": "bins must be an integer of at least 2 or a list of cut points, got {bins!r}.",
            "fr": "bins doit être un entier d'au moins 2 ou une liste de bornes, reçu {bins!r}.",
        },
        "describe.xr.bridged": {
            "en": "{count} of the {total} transitions bridge a gap in the calendar.",
            "fr": "{count} des {total} transitions enjambent un trou du calendrier.",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class XtTrans(Result):
    """Result of :func:`xttrans`.

    Attributes
    ----------
    variable
        Name of the variable.
    counts
        Number of transitions: rows are the value at the earlier observation, columns the value
        at the next one.
    percent
        The same in percent of each row.
    n_trans
        Transitions counted (``r(n_trans)``).
    n_rows, n_cols
        Distinct values of the rows and of the columns (``r(n_rows)``, ``r(n_cols)``).
    n_bridged
        Counted (``gaps="bridge"``) or discarded (``gaps="strict"``) transitions that span a gap
        in the calendar.
    gaps
        ``"bridge"`` (Stata) or ``"strict"``.
    """

    variable: str
    counts: pd.DataFrame
    percent: pd.DataFrame
    n_trans: int
    n_rows: int
    n_cols: int
    n_bridged: int
    gaps: Literal["bridge", "strict"] = "bridge"

    def to_frame(self) -> pd.DataFrame:
        """Return the row percentages."""
        return self.percent

    def summary(self, lang: Lang | None = None) -> str:
        """Return the matrix in the layout of Stata's ``xttrans`` (percent of each row)."""
        header = " " * 10 + "|" + "".join(f"{number(float(c), 6):>9}" for c in self.percent.columns)
        lines = [f"{self.variable:>10} |  {self.variable}", header + "    Total"]
        lines.append("-" * 11 + "+" + "-" * (9 * len(self.percent.columns) + 9))
        for key, row in self.percent.iterrows():
            cells = "".join(f"{v:>9.2f}" for v in row)
            lines.append(f"{number(float(str(key)), 6):>10} |{cells}   100.00")
        total = self.counts.to_numpy().sum(axis=0) / max(self.n_trans, 1) * 100
        lines.append("-" * 11 + "+" + "-" * (9 * len(self.percent.columns) + 9))
        lines.append(f"{'Total':>10} |" + "".join(f"{v:>9.2f}" for v in total) + "   100.00")
        if self.n_bridged:
            lines.append(
                _i18n.message("describe.xr.bridged", lang, count=self.n_bridged, total=self.n_trans)
            )
        return "\n".join(lines)


def discretize(x: FloatArray, bins: int | Sequence[float]) -> FloatArray:
    """Class number (``0`` lowest) of each value; classes are ``[c_j, c_{j+1})``, last one closed.

    An integer ``bins`` cuts at the quantiles ``1/bins, 2/bins, ...`` of the observed values;
    equal cut points (ties) merge classes.
    """
    observed = x[~np.isnan(x)]
    cuts: FloatArray
    if isinstance(bins, (int, np.integer)):
        if bins < 2:
            raise SpecificationError(
                _i18n.message("describe.xr.bad_bins", bins=bins), key="describe.xr.bad_bins"
            )
        cuts = (
            np.unique(np.quantile(observed, np.arange(1, bins) / bins))
            if observed.size
            else np.empty(0)
        )
    else:
        cuts = np.unique(np.asarray(list(bins), dtype=float))
    out = np.searchsorted(cuts, x, side="right").astype(np.float64)
    out[np.isnan(x)] = np.nan
    return out


@description(
    "xttrans",
    DescribeSpec(
        stata="xttrans",
        summary="Transition matrix between a unit's consecutive observations.",
    ),
)
def xttrans(
    panel: PanelBase,
    variable: str,
    *,
    gaps: Literal["bridge", "strict"] = "bridge",
    bins: int | Sequence[float] | None = None,
    max_values: int = 1000,
) -> XtTrans:
    """Count the transitions of a variable between consecutive observations of each unit.

    Parameters
    ----------
    panel
        The declared panel.
    variable
        A column (or a Stata-style term).  Like Stata, values are compared as ``float32``.
    gaps
        ``"bridge"`` (default, Stata): pair each observation with the unit's next row, whatever
        the calendar gap.  ``"strict"``: only pairs of consecutive periods.
    bins
        Discretize a continuous variable first (an addition of the package, not in Stata): an
        integer ``k`` cuts it at its ``k``-quantiles, a list gives the cut points.  Values are
        replaced by the number of the class, ``0`` for the lowest.
    max_values
        Refuse variables with more distinct values than this.

    Returns
    -------
    XtTrans
        Result of :func:`xttrans`.

    Raises
    ------
    DataError
        If the variable has no transition.
    SpecificationError
        If ``gaps`` is invalid or the variable has too many distinct values.
    """
    if gaps not in ("bridge", "strict"):
        raise SpecificationError(
            _i18n.message("describe.xr.bad_gaps", gaps=gaps), key="describe.xr.bad_gaps"
        )
    names, values = term_columns(panel, variable)
    if len(names) != 1:
        raise SpecificationError(
            _i18n.message("describe.xt.one_variable", command="xttrans", count=len(names)),
            key="describe.xt.one_variable",
        )
    # STATA-QUIRK: xttrans pairs a row with the unit's next row, so calendar gaps are bridged
    # STATA-QUIRK: xttrans stores the variable as float (single precision) before tabulating
    x = values[:, 0] if bins is None else discretize(values[:, 0], bins)
    x = x.astype(np.float32).astype(np.float64)
    index = panel.index
    first, second = np.arange(index.n_obs - 1), np.arange(1, index.n_obs)
    pair = (index.unit[first] == index.unit[second]) & ~np.isnan(x[first]) & ~np.isnan(x[second])
    spans_gap = index.rank[second] - index.rank[first] > 1
    n_bridged = int(np.count_nonzero(pair & spans_gap))
    if gaps == "strict":
        pair &= ~spans_gap
    if not pair.any():
        raise DataError(
            _i18n.message("describe.xr.no_transition", name=names[0]),
            stata_code=2000,
            key="describe.xr.no_transition",
        )
    was, now = x[first][pair], x[second][pair]
    rows, row_code = np.unique(was, return_inverse=True)
    cols, col_code = np.unique(now, return_inverse=True)
    if max(rows.size, cols.size) > max_values:
        raise SpecificationError(
            _i18n.message(
                "describe.xt.too_many",
                name=names[0],
                count=max(rows.size, cols.size),
                limit=max_values,
            ),
            stata_code=134,
            key="describe.xt.too_many",
        )
    table = np.bincount(row_code * cols.size + col_code, minlength=rows.size * cols.size).reshape(
        rows.size, cols.size
    )
    counts = pd.DataFrame(
        table, index=pd.Index(rows, name="from"), columns=pd.Index(cols, name="to")
    )
    percent = counts.div(counts.sum(axis=1), axis=0) * 100.0
    return XtTrans(
        variable=names[0],
        counts=counts,
        percent=percent,
        n_trans=int(table.sum()),
        n_rows=int(rows.size),
        n_cols=int(cols.size),
        n_bridged=n_bridged,
        gaps=gaps,
    )
