"""Atypical values inside units and level shifts, marked and never removed (K.4).

Additions of the package (the closest tools flag outliers across the whole sample, not within
each unit); validated by hand on small panels, against an independent loop implementation and by
invariance tests.

**Within outliers.**  For each unit, the *modified z-score* of Iglewicz and Hoaglin
(1993) compares a value with the median of the unit in units of the median absolute deviation::

    z_it = 0.6745 * (x_it - median_i) / MAD_i

where ``MAD_i`` is the median of ``|x_it - median_i|``.  When ``MAD_i`` is zero (more than half
of the unit's values are equal) the mean absolute deviation replaces it:
``z_it = (x_it - median_i) / (1.253314 * meanAD_i)``; a unit whose values are all equal gets
``z = 0``.  A unit with fewer than ``min_obs`` values is not scored.  Values with ``|z|`` above
``threshold`` (3.5 by default, the usual cut-off) are flagged.  The score is noisy when a unit
has few values: on normal data, simulated, 5.2 % of the values are flagged with 5 values per
unit, 3.2 % with 6, 2.4 % with 8, 1.3 % with 12, 0.7 % with 20 and 0.3 % with 50 (nominal 0.05 %),
so a small share of flagged values in short units is not evidence of outliers.

**Level shifts.**  In the sequence of the observed values of a unit, the shift before position
``j`` is the median of the next ``window`` values minus the median of the previous ``window``
values.  It is divided by the robust scale of the unit's first differences,
``1.4826 * MAD(diff)`` (the median of the positive scales of the other units when that is zero).
A spike moves one value and leaves both medians nearly unchanged; a lasting change of level
moves one of them.  The largest absolute score of each unit is reported.  A unit is
*influential* when its share of the total within sum of squares exceeds ``influence`` times the
average share.  Gaps between observations are ignored: the sequence is that of the observed values.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.explore._spec import ToolSpec, tool
from pyxtpanel.panel.arrays import sorted_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "explore",
    {
        "explore.out.title": {
            "en": "Within outliers of {variable}: {flagged} of {scored} scored values flagged "
            "(|z| > {threshold})",
            "fr": "Atypiques within de {variable} : {flagged} valeurs sur {scored} évaluées "
            "(|z| > {threshold})",
        },
        "explore.out.none": {"en": "No value flagged.", "fr": "Aucune valeur signalée."},
        "explore.shift.title": {
            "en": "Level shifts of {variable} (window {window}): {flagged} of {scored} units "
            "above {threshold}; {influential} influential units",
            "fr": "Sauts de niveau de {variable} (fenêtre {window}) : {flagged} unités sur "
            "{scored} au-dessus de {threshold} ; {influential} unités influentes",
        },
        "explore.out.args": {
            "en": "{name} must be {rule}, got {value!r}.",
            "fr": "{name} doit être {rule}, reçu {value!r}.",
        },
    },
)

#: Constants of the modified z-score.
_Z_SCALE, _MEAN_AD_SCALE, _MAD_TO_SD = 0.6745, 1.253314, 1.4826


def _check(name: str, ok: bool, rule: str, value: object) -> None:
    if not ok:
        raise SpecificationError(
            _i18n.message("explore.out.args", name=name, rule=rule, value=value),
            stata_code=198,
            key="explore.out.args",
        )


@dataclass(frozen=True, kw_only=True)
class WithinOutliers(Result):
    """Result of :func:`within_outliers`.

    Attributes
    ----------
    scores
        Modified z-score of each row, on the index of the declared frame (missing where the
        value is missing or the unit has fewer than ``min_obs`` values).
    flagged
        ``True`` where ``|score| > threshold``, on the same index.
    rows
        The flagged rows, largest ``|score|`` first: ``unit``, ``time``, ``value``, ``median``
        and ``score``.
    by_unit
        One row per unit: ``n_scored``, ``n_flagged``.
    variable, threshold, min_obs
        The specification.
    n_scored, n_flagged
        Totals.
    """

    scores: pd.Series
    flagged: pd.Series
    rows: pd.DataFrame
    by_unit: pd.DataFrame
    variable: str
    threshold: float
    min_obs: int
    n_scored: int
    n_flagged: int

    def to_frame(self) -> pd.DataFrame:
        """Return the flagged rows."""
        return self.rows

    def summary(self, lang: Lang | None = None) -> str:
        """Return the totals and the flagged rows as text."""
        head = _i18n.message(
            "explore.out.title",
            lang,
            variable=self.variable,
            flagged=self.n_flagged,
            scored=self.n_scored,
            threshold=self.threshold,
        )
        if not len(self.rows):
            return head + "\n" + _i18n.message("explore.out.none", lang)
        return head + "\n\n" + self.rows.head(30).to_string(float_format=lambda v: f"{v:.3f}")


@tool(
    "within_outliers",
    ToolSpec(
        summary="Within-unit outliers by robust scores (median and MAD in the unit), flagged "
        "and never removed.",
        validation="hand-computed panels, independent loop implementation, invariance tests",
    ),
)
def within_outliers(
    panel: PanelBase, variable: str, *, threshold: float = 3.5, min_obs: int = 5
) -> WithinOutliers:
    """Flag the values that are far from the median of their own unit.

    Parameters
    ----------
    panel
        A declared panel.
    variable
        The variable to examine.
    threshold
        A value is flagged when the absolute modified z-score exceeds it.
    min_obs
        Units with fewer observed values are not scored.
    """
    _check("threshold", threshold > 0, "positive", threshold)
    _check(
        "min_obs", isinstance(min_obs, int) and min_obs >= 2, "an integer of at least 2", min_obs
    )
    index = panel.index
    x = sorted_columns(panel, [variable])[:, 0]
    # Stata: egen median() averages the two middle values in float precision (1e-8 off); the exact
    # median used here agrees with summarize, detail and egen pctile(, p(50))
    series = pd.Series(x)
    by = series.groupby(index.unit)
    count = by.transform("count").to_numpy()
    median = by.transform("median").to_numpy()
    deviation = (series - median).abs()
    mad = deviation.groupby(index.unit).transform("median").to_numpy()
    mean_ad = deviation.groupby(index.unit).transform("mean").to_numpy()
    centred = x - median
    with np.errstate(invalid="ignore", divide="ignore"):
        score = np.where(
            mad > 0,
            _Z_SCALE * centred / mad,
            np.where(mean_ad > 0, centred / (_MEAN_AD_SCALE * mean_ad), 0.0),
        )
    score = np.where(np.isnan(x) | (count < min_obs), np.nan, score)
    scored = ~np.isnan(score)
    flag = scored & (np.abs(np.where(scored, score, 0.0)) > threshold)
    order = np.flatnonzero(flag)
    order = order[np.argsort(-np.abs(score[order]), kind="stable")]
    rows = pd.DataFrame(
        {
            "unit": index.unit_labels[index.unit[order]],
            "time": panel.axis.values(index.rank[order]),
            "value": x[order],
            "median": median[order],
            "score": score[order],
        }
    )
    by_unit = pd.DataFrame(
        {
            "n_scored": np.bincount(index.unit[scored], minlength=index.n_units),
            "n_flagged": np.bincount(index.unit[flag], minlength=index.n_units),
        },
        index=pd.Index(index.unit_labels, name=panel.id_name),
    )
    return WithinOutliers(
        scores=cast("pd.Series", panel.align(score, name="z")),
        flagged=cast("pd.Series", panel.align(flag, name="flagged")).fillna(False).astype(bool),
        rows=rows,
        by_unit=by_unit,
        variable=variable,
        threshold=float(threshold),
        min_obs=int(min_obs),
        n_scored=int(scored.sum()),
        n_flagged=int(flag.sum()),
    )


@dataclass(frozen=True, kw_only=True)
class LevelShifts(Result):
    """Result of :func:`level_shifts`.

    Attributes
    ----------
    by_unit
        One row per unit: ``n_obs``, ``max_score`` (absolute), ``shift`` (signed change of
        level at the largest score), ``period`` (first period after the shift), ``flagged``,
        ``within_share`` (share of the total within sum of squares) and ``influential``.
    variable, window, threshold, influence
        The specification.
    n_flagged, n_influential, n_scored
        Units flagged, influential and scored.
    """

    by_unit: pd.DataFrame
    variable: str
    window: int
    threshold: float
    influence: float
    n_flagged: int
    n_influential: int
    n_scored: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table by unit."""
        return self.by_unit

    def summary(self, lang: Lang | None = None) -> str:
        """Return the totals and the units flagged or influential as text."""
        head = _i18n.message(
            "explore.shift.title",
            lang,
            variable=self.variable,
            window=self.window,
            flagged=self.n_flagged,
            scored=self.n_scored,
            threshold=self.threshold,
            influential=self.n_influential,
        )
        table = self.by_unit
        shown = table[table["flagged"] | table["influential"]]
        if not len(shown):
            return head
        shown = shown.sort_values("max_score", ascending=False).head(30)
        return head + "\n\n" + shown.to_string(float_format=lambda v: f"{v:.3f}")


@tool(
    "level_shifts",
    ToolSpec(
        summary="Level shifts of each unit (change of median) and units that weigh most in the "
        "within variation.",
        validation="hand-computed panels, independent loop implementation, invariance tests",
    ),
)
def level_shifts(
    panel: PanelBase,
    variable: str,
    *,
    window: int = 3,
    threshold: float = 5.0,
    influence: float = 3.0,
) -> LevelShifts:
    """Find the largest change of level of each unit and the units that dominate the variation.

    Parameters
    ----------
    panel
        A declared panel.
    variable
        The variable to examine.
    window
        Number of observed values on each side of a candidate shift.
    threshold
        A unit is flagged when its largest absolute score exceeds it.
    influence
        A unit is influential when its share of the within sum of squares is above this
        multiple of the average share.
    """
    _check("window", isinstance(window, int) and window >= 1, "an integer of at least 1", window)
    _check("threshold", threshold > 0, "positive", threshold)
    _check("influence", influence > 0, "positive", influence)
    index = panel.index
    x = sorted_columns(panel, [variable])[:, 0]
    n_units = index.n_units
    n_obs = np.zeros(n_units, dtype=np.int64)
    scale = np.full(n_units, np.nan)
    best_shift = np.full(n_units, np.nan)
    best_rank = np.full(n_units, -1, dtype=np.int64)
    ss = np.zeros(n_units)
    candidates: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for unit in range(n_units):
        lo, hi = index.offsets[unit], index.offsets[unit + 1]
        keep = ~np.isnan(x[lo:hi])
        v, ranks = x[lo:hi][keep], index.rank[lo:hi][keep]
        n_obs[unit] = v.size
        if v.size:
            ss[unit] = float(np.sum((v - v.mean()) ** 2))
        if v.size < 2 * window or v.size < 2:
            continue
        jumps = np.diff(v)
        spread = _MAD_TO_SD * float(np.median(np.abs(jumps - np.median(jumps))))
        scale[unit] = spread
        medians = np.median(sliding_window_view(v, window), axis=1)
        shift = medians[window:] - medians[: v.size - 2 * window + 1]
        candidates[unit] = (shift, ranks[window : v.size - window + 1])
    positive = scale[scale > 0]
    fallback = float(np.median(positive)) if positive.size else np.nan
    score = np.full(n_units, np.nan)
    for unit, (shift, ranks) in candidates.items():
        used = scale[unit] if scale[unit] > 0 else fallback
        if not used > 0:
            continue
        # candidate positions are those with a full window on each side
        ratio = shift / used
        pick = int(np.argmax(np.abs(ratio)))
        score[unit] = abs(ratio[pick])
        best_shift[unit] = shift[pick]
        best_rank[unit] = ranks[pick]
    scored = ~np.isnan(score)
    total = ss.sum()
    share = ss / total if total > 0 else np.zeros(n_units)
    eligible = n_obs >= 2
    average = 1.0 / max(int(eligible.sum()), 1)
    influential = eligible & (share > influence * average)
    flagged = scored & (score > threshold)
    by_unit = pd.DataFrame(
        {
            "n_obs": n_obs,
            "max_score": score,
            "shift": best_shift,
            "period": pd.Series(panel.axis.values(np.maximum(best_rank, 0)), index=range(n_units))
            .where(best_rank >= 0)
            .to_numpy(),
            "flagged": flagged,
            "within_share": share,
            "influential": influential,
        },
        index=pd.Index(index.unit_labels, name=panel.id_name),
    )
    return LevelShifts(
        by_unit=by_unit,
        variable=variable,
        window=int(window),
        threshold=float(threshold),
        influence=float(influence),
        n_flagged=int(flagged.sum()),
        n_influential=int(influential.sum()),
        n_scored=int(scored.sum()),
    )
