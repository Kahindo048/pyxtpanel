"""Thresholds of the audit and the shared context its checks read (layer 4)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields, replace
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel.describe.attrition import Attrition, attrition
from pyxtpanel.describe.missing import MissingReport, missing_report
from pyxtpanel.explore.balance import Unbalancedness, unbalancedness
from pyxtpanel.explore.equations import Transform, UsableEquations, usable_equations
from pyxtpanel.explore.gaps import GapTypes, gap_types
from pyxtpanel.explore.outliers import WithinOutliers, within_outliers
from pyxtpanel.explore.variation import residual_variation, variation
from pyxtpanel.panel.arrays import sorted_columns

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "explore",
    {
        "explore.audit.threshold": {
            "en": "unknown audit threshold {name!r}; available: {available}.",
            "fr": "seuil d'audit inconnu {name!r} ; disponibles : {available}.",
        },
    },
)


@dataclass(frozen=True)
class AuditThresholds:
    """Limits above or below which the checks raise an alert.

    Override some of them with ``audit(panel, thresholds={"min_units": 50})``.

    Attributes
    ----------
    min_units
        Fewer units than this: ``few_units``.
    min_periods
        A calendar shorter than this: ``short_panel``.
    single_obs_share
        Share of units with one observation above which that is a warning.
    gamma
        Ahrens-Pincus gamma below this: ``unbalanced``.
    gap_share
        Share of units with an internal gap above which that is a warning.
    exit_share
        Share of units that leave early (or enter late) above which it is flagged.
    missing_share, missing_critical
        Share of missing values of a variable above which it is flagged / serious.
    complete_share
        Share of complete rows below which listwise deletion is a warning.
    within_share
        Share of the variance left by the unit-mean removal below which the variable has
        little within variation.
    outlier_threshold, outlier_share, outlier_min_obs
        Modified z-score cut-off, share of flagged values that makes it a warning, and
        number of values a unit needs to be scored.  Scores are noisy in short units: on normal
        data 3.2 % of the values are flagged at 3.5 with 6 per unit, 1.3 % with 12, 0.3 % with 50,
        hence the 5 % default.
    equations_share
        Share of the observations that give a first-difference equation below which it is a
        warning.
    fod_gain
        Relative gain of orthogonal deviations over first differences that makes them
        recommended.
    """

    min_units: int = 30
    min_periods: int = 4
    single_obs_share: float = 0.10
    gamma: float = 0.80
    gap_share: float = 0.10
    exit_share: float = 0.25
    missing_share: float = 0.20
    missing_critical: float = 0.50
    complete_share: float = 0.50
    within_share: float = 0.05
    outlier_threshold: float = 3.5
    outlier_share: float = 0.05
    outlier_min_obs: int = 6
    equations_share: float = 0.50
    fod_gain: float = 0.20

    @classmethod
    def from_mapping(cls, values: Mapping[str, float] | None) -> AuditThresholds:
        """Defaults with ``values`` applied; an unknown name raises ``SpecificationError``."""
        base = cls()
        if not values:
            return base
        known = {f.name for f in fields(cls)}
        for name in values:
            if name not in known:
                raise SpecificationError(
                    _i18n.message(
                        "explore.audit.threshold", name=name, available=", ".join(sorted(known))
                    ),
                    stata_code=198,
                    key="explore.audit.threshold",
                )
        return replace(base, **dict(values))  # type: ignore[arg-type]


class AuditContext:
    """A panel, its variables and the results the checks share (each computed once)."""

    def __init__(
        self,
        panel: PanelBase,
        variables: Sequence[str],
        thresholds: AuditThresholds,
        dependent: str | None = None,
        lags: int = 1,
    ) -> None:
        self.panel = panel
        self.variables = list(variables)
        self.thresholds = thresholds
        self.dependent = dependent
        self.lags = lags
        self._outliers: dict[str, WithinOutliers] = {}
        self._equations: dict[str, UsableEquations] = {}

    @cached_property
    def gaps(self) -> GapTypes:
        """Typology of the gaps."""
        return gap_types(self.panel)

    @cached_property
    def balance(self) -> Unbalancedness:
        """Imbalance measures."""
        return unbalancedness(self.panel)

    @cached_property
    def attrition(self) -> Attrition:
        """Entries and exits."""
        return attrition(self.panel)

    @cached_property
    def missing(self) -> MissingReport | None:
        """Missing values (``None`` without variables)."""
        return missing_report(self.panel, self.variables) if self.variables else None

    @cached_property
    def complete_share(self) -> float:
        """Share of the rows where every variable is observed."""
        if not self.variables or self.panel.n_obs == 0:
            return 1.0
        data = sorted_columns(self.panel, self.variables)
        return float((~np.isnan(data).any(axis=1)).mean())

    @cached_property
    def observed_variables(self) -> list[str]:
        """Variables with at least one observed value (the others are reported as missing)."""
        if not self.variables:
            return []
        data = sorted_columns(self.panel, self.variables)
        return [
            v for v, ok in zip(self.variables, (~np.isnan(data)).any(axis=0), strict=True) if ok
        ]

    @cached_property
    def variation(self) -> pd.DataFrame:
        """Table of :func:`~pyxtpanel.explore.variation.variation`."""
        names = self.observed_variables
        return variation(self.panel, names).table if names else pd.DataFrame()

    @cached_property
    def within(self) -> pd.DataFrame:
        """Variance left by the within transformation, by variable."""
        names = self.observed_variables
        if not names:
            return pd.DataFrame()
        table = residual_variation(self.panel, names, ["within"]).table
        return table.droplevel("transform")

    def outliers(self, variable: str) -> WithinOutliers:
        """Within outliers of ``variable`` (cached)."""
        if variable not in self._outliers:
            self._outliers[variable] = within_outliers(
                self.panel,
                variable,
                threshold=self.thresholds.outlier_threshold,
                min_obs=self.thresholds.outlier_min_obs,
            )
        return self._outliers[variable]

    def equations(self, transform: Transform) -> UsableEquations:
        """Usable equations of the dependent variable (cached by transformation)."""
        if self.dependent is None:
            raise ValueError("no dependent variable was given")
        if transform not in self._equations:
            self._equations[transform] = usable_equations(
                self.panel, self.dependent, lags=self.lags, transform=transform
            )
        return self._equations[transform]
