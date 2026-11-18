"""Measures of the imbalance of a panel (K.2.2).

Same two measures as ``punbalancedness`` of the R package ``plm`` (Ahrens and Pincus, 1981), for
a panel with units ``i = 1..N`` observed ``T_i`` times, ``n = sum(T_i)`` observations and
``Tbar = n / N``:

* ``gamma = N / (Tbar * sum(1 / T_i))``, the harmonic mean of the ``T_i`` over their mean;
* ``nu = 1 / (N * sum((T_i / n)^2))``, the inverse of ``N`` times the sum of the squared shares
  of the observations.

Both equal 1 for a balanced panel and fall towards 0 as the panel becomes less balanced
(``1 / N < nu <= 1``).  The result also gives the plain figures behind them (``Tbar``,
harmonic mean, fill rate of the calendar).  Validated against the published formulas, by hand
on small panels and against an independent loop implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.explore._spec import ToolSpec, tool

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "explore",
    {
        "explore.balance.title": {
            "en": "Imbalance of the panel: {n_units} units, {n_obs} observations",
            "fr": "Déséquilibre du panel : {n_units} unités, {n_obs} observations",
        },
        "explore.balance.balanced": {
            "en": "The panel is balanced: every unit is observed in all {n} periods.",
            "fr": "Le panel est cylindré : chaque unité est observée aux {n} périodes.",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class Unbalancedness(Result):
    """Result of :func:`unbalancedness`.

    Attributes
    ----------
    gamma, nu
        The two measures of Ahrens and Pincus; 1 when balanced.
    n_units, n_obs, n_periods
        Units, observations and periods of the calendar of the panel.
    t_mean, t_harmonic, t_min, t_max
        Mean, harmonic mean, minimum and maximum of ``T_i``.
    fill_rate
        ``n_obs / (n_units * n_periods)``: share of the cells of the full grid observed.
    """

    gamma: float
    nu: float
    n_units: int
    n_obs: int
    n_periods: int
    t_mean: float
    t_harmonic: float
    t_min: int
    t_max: int
    fill_rate: float

    def to_frame(self) -> pd.DataFrame:
        """Return the measures as a one-column table."""
        values = {
            "gamma": self.gamma,
            "nu": self.nu,
            "T_mean": self.t_mean,
            "T_harmonic": self.t_harmonic,
            "T_min": self.t_min,
            "T_max": self.t_max,
            "fill_rate": self.fill_rate,
        }
        # a Series keeps the order of the keys on every pandas version (a dict of dicts is sorted
        # by pandas before 2.2)
        return pd.Series(values, name="value").to_frame()

    def summary(self, lang: Lang | None = None) -> str:
        """Return the measures as text."""
        lines = [
            _i18n.message("explore.balance.title", lang, n_units=self.n_units, n_obs=self.n_obs),
            "",
            self.to_frame().to_string(float_format=lambda v: f"{v:.4f}"),
        ]
        if self.gamma == 1.0 and self.nu == 1.0:
            lines += ["", _i18n.message("explore.balance.balanced", lang, n=self.n_periods)]
        return "\n".join(lines)


@tool(
    "unbalancedness",
    ToolSpec(
        summary="Ahrens-Pincus measures of the imbalance of a panel (gamma and nu).",
        validation="published formulas, hand-computed panels, independent loop implementation",
        similar="plm::punbalancedness",
    ),
)
def unbalancedness(panel: PanelBase) -> Unbalancedness:
    """Compute ``gamma`` and ``nu`` (as ``plm::punbalancedness``) and the figures behind them."""
    counts = panel.index.counts.astype(np.float64)
    n_units, n_obs = counts.size, float(counts.sum())
    t_mean = n_obs / n_units
    gamma = n_units / (t_mean * float(np.sum(1.0 / counts)))
    nu = 1.0 / (n_units * float(np.sum((counts / n_obs) ** 2)))
    n_periods = int(np.max(panel.index.rank) - np.min(panel.index.rank) + 1)
    return Unbalancedness(
        gamma=gamma,
        nu=nu,
        n_units=n_units,
        n_obs=int(n_obs),
        n_periods=n_periods,
        t_mean=t_mean,
        t_harmonic=n_units / float(np.sum(1.0 / counts)),
        t_min=int(np.min(counts)),
        t_max=int(np.max(counts)),
        fill_rate=n_obs / (n_units * n_periods),
    )


@dataclass(frozen=True)
class Block:
    """A balanced block: the units observed in every period of a window of consecutive periods.

    Attributes
    ----------
    start, end
        First and last period of the window (the user's time values).
    n_periods, n_units, n_obs
        Length of the window, units complete in it and observations (``n_units * n_periods``).
    """

    start: object
    end: object
    n_periods: int
    n_units: int
    n_obs: int


@dataclass(frozen=True, kw_only=True)
class BalanceCost(Result):
    """Result of :func:`balance_cost`.

    Attributes
    ----------
    thresholds
        Indexed by ``k`` (``T_i >= k``): ``units`` and ``obs`` kept (all the observations of
        those units) and their shares of the panel.  The panel is not balanced by this rule; it
        shows what a minimum length costs.
    windows
        Indexed by ``length`` of a window of consecutive periods: the earliest ``start`` that
        keeps the most units, ``end``, ``units`` complete in it, ``obs`` (``units * length``) and
        the shares of the panel.  ``units`` never increases with ``length``.
    largest
        The window and units that hold the most observations (what
        ``PanelData.balance("largest_block")`` keeps).
    n_units, n_obs, n_periods
        Size of the panel and of its calendar.
    """

    thresholds: pd.DataFrame
    windows: pd.DataFrame
    largest: Block
    n_units: int
    n_obs: int
    n_periods: int

    def to_frame(self) -> pd.DataFrame:
        """Return the table by window length."""
        return self.windows

    def summary(self, lang: Lang | None = None) -> str:
        """Return both tables and the largest block as text."""
        big = self.largest
        return "\n".join(
            [
                _i18n.message(
                    "explore.cost.title",
                    lang,
                    n_units=self.n_units,
                    n_obs=self.n_obs,
                    n_periods=self.n_periods,
                ),
                "",
                _i18n.message("explore.cost.thresholds", lang),
                self.thresholds.to_string(float_format=lambda v: f"{v:.3f}"),
                "",
                _i18n.message("explore.cost.windows", lang),
                self.windows.to_string(float_format=lambda v: f"{v:.3f}"),
                "",
                _i18n.message(
                    "explore.cost.largest",
                    lang,
                    start=big.start,
                    end=big.end,
                    units=big.n_units,
                    periods=big.n_periods,
                    obs=big.n_obs,
                    share=100.0 * big.n_obs / self.n_obs,
                ),
            ]
        )


_i18n.register(
    "explore",
    {
        "explore.cost.title": {
            "en": "Cost of balancing: {n_units} units, {n_obs} observations, {n_periods} periods",
            "fr": "Coût du cylindrage : {n_units} unités, {n_obs} observations, "
            "{n_periods} périodes",
        },
        "explore.cost.thresholds": {
            "en": "Units kept if T_i >= k (all their observations):",
            "fr": "Unités conservées si T_i >= k (toutes leurs observations) :",
        },
        "explore.cost.windows": {
            "en": "Best window of consecutive periods, by length (units complete in it):",
            "fr": "Meilleure fenêtre de périodes consécutives, par longueur (unités complètes) :",
        },
        "explore.cost.largest": {
            "en": "Largest balanced block: {start} to {end}, {units} units x {periods} periods = "
            "{obs} observations ({share:.1f} % of the panel).",
            "fr": "Plus grand bloc cylindré : {start} à {end}, {units} unités x {periods} périodes "
            "= {obs} observations ({share:.1f} % du panel).",
        },
    },
)


@tool(
    "balance_cost",
    ToolSpec(
        summary="Units and observations kept by each minimum length T_i >= k and each window.",
        validation="hand-computed panels, independent loop implementation, "
        "same block as PanelData.balance('largest_block')",
        similar="Stata xtbalfront, xtbalance2 (SSC); plm::make.pbalanced applies one rule",
    ),
)
def balance_cost(panel: PanelBase) -> BalanceCost:
    """Show what each way of balancing the panel keeps, before choosing one.

    Two families of rules are tabulated: keeping the units with at least ``k`` observations (the
    result is still unbalanced), and keeping a window of consecutive periods with the units
    complete in it (the result is balanced).  For each window length the best start is the one
    that keeps the most units, the earliest in case of ties.  Only windows of consecutive
    periods are searched: the largest rectangle with arbitrary periods is an NP-hard problem.
    """
    index = panel.index
    n_units, n_obs = index.n_units, index.n_obs
    low, high = int(np.min(index.rank)), int(np.max(index.rank))
    span = high - low + 1

    top = int(np.max(index.counts))
    ks = np.arange(1, top + 1)
    units_k = np.array([(index.counts >= k).sum() for k in ks])
    obs_k = np.array([index.counts[index.counts >= k].sum() for k in ks])
    thresholds = pd.DataFrame(
        {
            "units": units_k,
            "obs": obs_k,
            "units_share": units_k / n_units,
            "obs_share": obs_k / n_obs,
        },
        index=pd.Index(ks, name="k"),
    )

    # run[i, s]: consecutive periods observed from s onwards; unit i is complete in the window
    # [s, s + L - 1] exactly when run[i, s] >= L
    run = np.zeros((n_units, span + 1), dtype=np.int64)
    present = np.zeros((n_units, span), dtype=np.bool_)
    present[index.unit, index.rank - low] = True
    for s in range(span - 1, -1, -1):
        run[:, s] = np.where(present[:, s], run[:, s + 1] + 1, 0)
    best_units = np.zeros(span, dtype=np.int64)
    best_start = np.zeros(span, dtype=np.int64)
    for length in range(1, span + 1):
        per_start = (run[:, : span - length + 1] >= length).sum(axis=0)
        start = int(np.argmax(per_start))  # first maximum: the earliest start
        best_units[length - 1], best_start[length - 1] = per_start[start], start
    lengths = np.arange(1, span + 1)
    obs_w = best_units * lengths
    windows = pd.DataFrame(
        {
            "start": panel.axis.values(best_start + low),
            "end": panel.axis.values(best_start + low + lengths - 1),
            "units": best_units,
            "obs": obs_w,
            "units_share": best_units / n_units,
            "obs_share": obs_w / n_obs,
        },
        index=pd.Index(lengths, name="length"),
    )
    # the same tie-break as balance("largest_block"): most observations, longer, earlier
    score = [
        (int(o), int(length), -int(s))
        for o, length, s in zip(obs_w, lengths, best_start, strict=True)
    ]
    pick = max(range(span), key=score.__getitem__)
    largest = Block(
        start=windows["start"].iloc[pick],
        end=windows["end"].iloc[pick],
        n_periods=int(lengths[pick]),
        n_units=int(best_units[pick]),
        n_obs=int(obs_w[pick]),
    )
    return BalanceCost(
        thresholds=thresholds,
        windows=windows,
        largest=largest,
        n_units=n_units,
        n_obs=n_obs,
        n_periods=span,
    )
