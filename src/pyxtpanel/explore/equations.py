"""Equations a dynamic-panel estimator can use, by unit and period (K.2.4).

An addition of the package (no software reports it before estimation); validated by hand on
small panels, against an independent loop implementation and against the number of observations
of ``pyxtabond2``.

A GMM estimator for ``y_it = sum_l phi_l y_{i,t-l} + x_it b + u_i + e_it`` transforms the model
to remove ``u_i`` and keeps the rows where every transformed variable exists:

* **first differences** (``"fd"``): ``D.y_t``, ``D.L.y_t`` ... ``D.L{p}.y_t`` and ``D.x_t``
  need ``y`` at ``t, t-1, ..., t-p-1`` and ``x`` at ``t, t-1``;
* **forward orthogonal deviations** (``"fod"``): ``y*_t`` and the transformed lags need ``y_t``,
  ``y_{t-1}, ..., y_{t-p}`` and at least one later observation of each variable of the unit,
  whatever the gaps (the definition of ``xtabond2``).

The transformations are those of :mod:`pyxtpanel.ops`, so the rows counted here are the rows the
estimator would receive.  A row is *instrumented* when, in addition, ``y`` is observed at least
``instrument_lag`` periods earlier (2 for first differences, 1 for orthogonal deviations): the
first level that can serve as a GMM instrument.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import BoolArray, FloatArray, Lang
from pyxtpanel.explore._spec import ToolSpec, tool
from pyxtpanel.ops.diff import difference, fod
from pyxtpanel.ops.lags import lag
from pyxtpanel.panel.arrays import sorted_columns
from pyxtpanel.panel.index import PanelIndex

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

Transform = Literal["fd", "fod"]

_i18n.register(
    "explore",
    {
        "explore.eq.title": {
            "en": "Usable equations for {variable} with {lags} lag(s), {transform}: "
            "{usable} of {observed} observations, {units} of {n_units} units",
            "fr": "Équations utilisables pour {variable} avec {lags} retard(s), {transform} : "
            "{usable} sur {observed} observations, {units} unités sur {n_units}",
        },
        "explore.eq.instrumented": {
            "en": "{n} of them have a GMM instrument (y observed {m} or more periods earlier).",
            "fr": "{n} d'entre elles ont un instrument GMM (y observé {m} périodes avant ou plus).",
        },
        "explore.eq.compare": {
            "en": "Equations by number of lags and transformation:",
            "fr": "Équations selon le nombre de retards et la transformation :",
        },
        "explore.eq.transform": {
            "en": "transform must be 'fd' or 'fod', got {value!r}.",
            "fr": "transform doit valoir 'fd' ou 'fod', reçu {value!r}.",
        },
        "explore.eq.lags": {
            "en": "lags must be a non-negative integer, got {value!r}.",
            "fr": "lags doit être un entier positif ou nul, reçu {value!r}.",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class UsableEquations(Result):
    """Result of :func:`usable_equations`.

    Attributes
    ----------
    by_row
        One row per observation, indexed by (unit, period): ``observed`` (the variable is
        observed), ``usable`` (every transformed variable exists) and ``instrumented``.
    by_unit
        One row per unit: ``n_obs``, ``n_usable`` and ``n_instrumented``.
    by_period
        One row per calendar period: ``observed``, ``usable``, ``instrumented`` counts.
    compare
        Equations (``equations``) and units (``units``) for each transformation and number of
        lags, with the same regressors, to choose the specification.
    variable, lags, transform, exog, instrument_lag
        The specification.
    n_observed, n_usable, n_instrumented, n_units_usable, n_units
        Totals.
    """

    by_row: pd.DataFrame
    by_unit: pd.DataFrame
    by_period: pd.DataFrame
    compare: pd.DataFrame
    variable: str
    lags: int
    transform: Transform
    exog: tuple[str, ...]
    instrument_lag: int
    n_observed: int
    n_usable: int
    n_instrumented: int
    n_units_usable: int
    n_units: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table by period."""
        return self.by_period

    def summary(self, lang: Lang | None = None) -> str:
        """Return the totals, the table by period and the comparison as text."""
        return "\n".join(
            [
                _i18n.message(
                    "explore.eq.title",
                    lang,
                    variable=self.variable,
                    lags=self.lags,
                    transform=self.transform,
                    usable=self.n_usable,
                    observed=self.n_observed,
                    units=self.n_units_usable,
                    n_units=self.n_units,
                ),
                _i18n.message(
                    "explore.eq.instrumented",
                    lang,
                    n=self.n_instrumented,
                    m=self.instrument_lag,
                ),
                "",
                self.by_period.to_string(),
                "",
                _i18n.message("explore.eq.compare", lang),
                self.compare.to_string(),
            ]
        )


def _mask(
    index: PanelIndex, y: FloatArray, exog: FloatArray, lags: int, transform: Transform
) -> BoolArray:
    """Rows where ``y``, its ``lags`` lags and ``exog`` all exist after the transformation."""
    columns = [y] + [lag(index, y, k) for k in range(1, lags + 1)]
    columns += [exog[:, j] for j in range(exog.shape[1])]
    ok = np.ones(index.n_obs, dtype=np.bool_)
    for column in columns:
        moved = difference(index, column) if transform == "fd" else fod(index, column, scale=False)
        ok &= ~np.isnan(moved)
    return ok


def _instrumented(index: PanelIndex, y: FloatArray, usable: BoolArray, m: int) -> BoolArray:
    """Usable rows with ``y`` observed at least ``m`` periods earlier."""
    seen = ~np.isnan(y)
    earliest = np.full(index.n_units, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(earliest, index.unit[seen], index.rank[seen])
    reach = index.rank.astype(np.int64) - m
    return usable & (earliest[index.unit] <= reach)


@tool(
    "usable_equations",
    ToolSpec(
        summary="Equations left for a dynamic-panel estimator, by unit and period (FD or FOD).",
        validation="hand-computed panels, independent loop implementation, "
        "observations of pyxtabond2",
    ),
)
def usable_equations(
    panel: PanelBase,
    variable: str,
    *,
    lags: int = 1,
    transform: Transform = "fd",
    exog: Sequence[str] = (),
    instrument_lag: int | None = None,
) -> UsableEquations:
    """Count the rows a dynamic-panel estimator keeps for ``variable`` with ``lags`` lags.

    Parameters
    ----------
    panel
        A declared panel.
    variable
        The dependent variable.
    lags
        Number of lags of ``variable`` among the regressors (0 for a static model).
    transform
        ``"fd"`` (first differences) or ``"fod"`` (forward orthogonal deviations).
    exog
        Other regressors that must exist after the transformation.
    instrument_lag
        Oldest lag allowed to serve as an instrument; default 2 for ``"fd"`` and 1 for
        ``"fod"``.

    Returns
    -------
    UsableEquations
        Result of :func:`usable_equations`.

    Raises
    ------
    SpecificationError
        If ``transform`` is not ``"fd"`` or ``"fod"``, or ``lags`` is not a non-negative integer.
    """
    if transform not in ("fd", "fod"):
        raise SpecificationError(
            _i18n.message("explore.eq.transform", value=transform),
            stata_code=198,
            key="explore.eq.transform",
        )
    if isinstance(lags, bool) or not isinstance(lags, (int, np.integer)) or lags < 0:
        raise SpecificationError(
            _i18n.message("explore.eq.lags", value=lags), stata_code=198, key="explore.eq.lags"
        )
    m = (2 if transform == "fd" else 1) if instrument_lag is None else int(instrument_lag)
    index = panel.index
    names = [variable, *exog]
    data = sorted_columns(panel, names)
    y, extra = data[:, 0], data[:, 1:]

    usable = _mask(index, y, extra, lags, transform)
    instrumented = _instrumented(index, y, usable, m)
    observed = ~np.isnan(y)

    multi = pd.MultiIndex.from_arrays(
        [index.unit_labels[index.unit], panel.axis.values(index.rank)],
        names=[panel.id_name, panel.time_name],
    )
    by_row = pd.DataFrame(
        {"observed": observed, "usable": usable, "instrumented": instrumented}, index=multi
    )
    starts = index.offsets[:-1]
    by_unit = pd.DataFrame(
        {
            "n_obs": np.add.reduceat(observed.astype(np.int64), starts),
            "n_usable": np.add.reduceat(usable.astype(np.int64), starts),
            "n_instrumented": np.add.reduceat(instrumented.astype(np.int64), starts),
        },
        index=pd.Index(index.unit_labels, name=panel.id_name),
    )
    low, high = int(np.min(index.rank)), int(np.max(index.rank))
    size = high - low + 1

    def per_period(flag: BoolArray) -> np.ndarray:
        return np.bincount(index.rank[flag] - low, minlength=size)

    by_period = pd.DataFrame(
        {
            "observed": per_period(observed),
            "usable": per_period(usable),
            "instrumented": per_period(instrumented),
        },
        index=pd.Index(panel.axis.labels(np.arange(low, high + 1)), name=panel.time_name),
    )
    rows = []
    for kind in ("fd", "fod"):
        for p in range(0, max(3, lags) + 1):
            ok = _mask(index, y, extra, p, kind)
            rows.append((kind, p, int(ok.sum()), int(np.unique(index.unit[ok]).size)))
    compare = pd.DataFrame(rows, columns=["transform", "lags", "equations", "units"]).set_index(
        ["transform", "lags"]
    )
    return UsableEquations(
        by_row=by_row,
        by_unit=by_unit,
        by_period=by_period,
        compare=compare,
        variable=variable,
        lags=int(lags),
        transform=transform,
        exog=tuple(exog),
        instrument_lag=m,
        n_observed=int(observed.sum()),
        n_usable=int(usable.sum()),
        n_instrumented=int(instrumented.sum()),
        n_units_usable=int(np.unique(index.unit[usable]).size),
        n_units=index.n_units,
    )
