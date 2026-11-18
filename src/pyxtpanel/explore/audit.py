"""``audit()``: the whole panel examined in one call, with alerts and advice (K.1).

An addition of the package (the closest are ``explore_panel_structure`` of ``expdpy`` and
``diagnose()`` of ``pyxtdpthresh``); validated by hand on small panels, by the loop oracles of
the tools it aggregates, and by tests of the alerts on panels built to raise each of them.

The audit runs the registered checks of :mod:`pyxtpanel.explore.checks` on a shared context, so
no result is computed twice.  Each check returns :class:`~pyxtpanel.explore.checks._spec.Alert`
objects with a level (critical, warning, info), a message and an advice in French and English.
The report exports to text, Markdown, HTML (also shown by Jupyter) and LaTeX.
"""

from __future__ import annotations

import html
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._format import latex_escape, markdown_table
from pyxtpanel._results import Result
from pyxtpanel._types import Lang
from pyxtpanel.explore._spec import ToolSpec, tool
from pyxtpanel.explore.checks._context import AuditContext, AuditThresholds
from pyxtpanel.explore.checks._spec import CHECKS, LEVELS, Alert, Level, run
from pyxtpanel.explore.variation import numeric_names

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "explore",
    {
        "explore.audit.title": {"en": "Panel audit", "fr": "Audit du panel"},
        "explore.audit.counts": {
            "en": "{critical} critical, {warning} warning(s), {info} for information",
            "fr": "{critical} critique(s), {warning} avertissement(s), {info} information(s)",
        },
        "explore.audit.level.critical": {"en": "CRITICAL", "fr": "CRITIQUE"},
        "explore.audit.level.warning": {"en": "WARNING", "fr": "ATTENTION"},
        "explore.audit.level.info": {"en": "INFO", "fr": "INFO"},
        "explore.audit.clean": {
            "en": "Nothing to report with these thresholds.",
            "fr": "Rien à signaler avec ces seuils.",
        },
        "explore.audit.facts": {"en": "Key figures", "fr": "Chiffres clés"},
        "explore.audit.variables": {"en": "Variables", "fr": "Variables"},
        "explore.audit.alerts": {"en": "Alerts", "fr": "Alertes"},
        "explore.audit.columns.level": {"en": "level", "fr": "niveau"},
        "explore.audit.columns.message": {"en": "finding", "fr": "constat"},
        "explore.audit.columns.advice": {"en": "advice", "fr": "conseil"},
        "explore.audit.check": {
            "en": "unknown audit check {name!r}; available: {available}.",
            "fr": "contrôle d'audit inconnu {name!r} ; disponibles : {available}.",
        },
        "explore.audit.duplicates": {
            "en": "{pairs} (unit, period) pairs appear several times ({rows} rows): the panel "
            "cannot be declared.",
            "fr": "{pairs} paires (unité, période) apparaissent plusieurs fois ({rows} lignes) : "
            "le panel ne peut pas être déclaré.",
        },
        "explore.audit.duplicates.advice": {
            "en": "Merge them with resolve_duplicates(frame, id, time, how=...) or list them "
            "with duplicates(frame, id, time), then audit again.",
            "fr": "Fusionnez-les avec resolve_duplicates(frame, id, time, how=...) ou listez-les "
            "avec duplicates(frame, id, time), puis relancez l'audit.",
        },
    },
)

_COLORS = {"critical": "#b3261e", "warning": "#b26a00", "info": "#2a78d6"}


@dataclass(frozen=True, kw_only=True)
class AuditReport(Result):
    """Result of :func:`audit`.

    Attributes
    ----------
    alerts
        Findings, most serious first.
    facts
        Key figures of the panel (``value`` column).
    variables
        One row per variable: ``n_obs``, ``pct_missing``, ``kind``, ``within_variance_share``
        and ``outliers_flagged``.
    n_critical, n_warning, n_info
        Number of alerts of each level.
    """

    alerts: tuple[Alert, ...]
    facts: pd.DataFrame
    variables: pd.DataFrame
    n_critical: int
    n_warning: int
    n_info: int

    def has_alerts(self, level: Level | None = None) -> bool:
        """Whether some alert exists (of the given level, if any)."""
        if level is None:
            return bool(self.alerts)
        return any(a.level == level for a in self.alerts)

    def codes(self) -> list[str]:
        """Codes of the alerts, in report order."""
        return [a.code for a in self.alerts]

    def to_frame(self, lang: Lang | None = None) -> pd.DataFrame:
        """Return the alerts as a table: level, check, code, message and advice."""
        return pd.DataFrame(
            {
                "level": [a.level for a in self.alerts],
                "check": [a.check for a in self.alerts],
                "code": [a.code for a in self.alerts],
                "message": [a.message(lang) for a in self.alerts],
                "advice": [a.advice(lang) for a in self.alerts],
            }
        )

    def _shown_facts(self) -> pd.DataFrame:
        """Key figures as text: counts without decimals, ratios with three."""
        shown = self.facts["value"].map(
            lambda v: f"{v:.0f}" if float(v).is_integer() else f"{v:.3f}"
        )
        return pd.DataFrame({"value": shown})

    def _heading(self, lang: Lang | None) -> str:
        counts = _i18n.message(
            "explore.audit.counts",
            lang,
            critical=self.n_critical,
            warning=self.n_warning,
            info=self.n_info,
        )
        return f"{_i18n.message('explore.audit.title', lang)}: {counts}"

    def summary(self, lang: Lang | None = None) -> str:
        """Return the report as text."""
        lines = [self._heading(lang), ""]
        if not self.alerts:
            lines.append(_i18n.message("explore.audit.clean", lang))
        for alert in self.alerts:
            tag = _i18n.message(f"explore.audit.level.{alert.level}", lang)
            lines.append(f"[{tag}] {alert.message(lang)}")
            lines.append(f"    -> {alert.advice(lang)}")
        if len(self.facts):
            lines += [
                "",
                _i18n.message("explore.audit.facts", lang),
                self._shown_facts().to_string(),
            ]
        if len(self.variables):
            lines += [
                "",
                _i18n.message("explore.audit.variables", lang),
                self.variables.to_string(float_format=lambda v: f"{v:.3f}"),
            ]
        return "\n".join(lines)

    def to_markdown(self, digits: int = 4, *, lang: Lang | None = None) -> str:
        """Return the report as Markdown (heading, alert list, tables)."""
        lines = [f"## {self._heading(lang)}", ""]
        if not self.alerts:
            lines.append(_i18n.message("explore.audit.clean", lang))
        for alert in self.alerts:
            tag = _i18n.message(f"explore.audit.level.{alert.level}", lang)
            lines.append(f"- **{tag}** {alert.message(lang)}  ")
            lines.append(f"  *{alert.advice(lang)}*")
        if len(self.facts):
            lines += ["", f"### {_i18n.message('explore.audit.facts', lang)}", ""]
            lines.append(markdown_table(self._shown_facts(), digits))
        if len(self.variables):
            lines += ["", f"### {_i18n.message('explore.audit.variables', lang)}", ""]
            lines.append(markdown_table(self.variables, digits))
        return "\n".join(lines)

    def to_latex(self, digits: int = 4, *, lang: Lang | None = None) -> str:
        """Return the alerts as a LaTeX ``tabular`` with wrapped text columns."""
        head = " & ".join(
            latex_escape(_i18n.message(f"explore.audit.columns.{c}", lang))
            for c in ("level", "message", "advice")
        )
        rows = [
            " & ".join(
                latex_escape(text)
                for text in (
                    _i18n.message(f"explore.audit.level.{a.level}", lang),
                    a.message(lang),
                    a.advice(lang),
                )
            )
            + r" \\"
            for a in self.alerts
        ]
        return "\n".join(
            [
                r"\begin{tabular}{lp{6.5cm}p{6.5cm}}",
                r"\hline",
                head + r" \\",
                r"\hline",
                *rows,
                r"\hline",
                r"\end{tabular}",
            ]
        )

    def to_html(self, lang: Lang | None = None) -> str:
        """Return the report as an HTML fragment (what Jupyter displays)."""
        parts = [f"<div class='pyxtpanel-audit'><h3>{html.escape(self._heading(lang))}</h3>"]
        if not self.alerts:
            parts.append(f"<p>{html.escape(_i18n.message('explore.audit.clean', lang))}</p>")
        for alert in self.alerts:
            tag = html.escape(_i18n.message(f"explore.audit.level.{alert.level}", lang))
            parts.append(
                f"<p style='border-left:4px solid {_COLORS[alert.level]};padding-left:8px'>"
                f"<b>{tag}</b> {html.escape(alert.message(lang))}<br>"
                f"<i>{html.escape(alert.advice(lang))}</i></p>"
            )
        if len(self.facts):
            parts.append(f"<h4>{html.escape(_i18n.message('explore.audit.facts', lang))}</h4>")
            parts.append(self._shown_facts().to_html())
        if len(self.variables):
            parts.append(f"<h4>{html.escape(_i18n.message('explore.audit.variables', lang))}</h4>")
            parts.append(self.variables.to_html(float_format=lambda v: f"{v:.3f}"))
        parts.append("</div>")
        return "".join(parts)

    def _repr_html_(self) -> str:
        return self.to_html()


def duplicates_report(n_rows: int, n_pairs: int) -> AuditReport:
    """Report for raw data that cannot be declared because (unit, period) pairs repeat."""
    alert = Alert("duplicates", "duplicates", "critical", {"rows": n_rows, "pairs": n_pairs})
    return AuditReport(
        alerts=(alert,),
        facts=pd.DataFrame({"value": []}, index=pd.Index([], name="fact")),
        variables=pd.DataFrame(),
        n_critical=1,
        n_warning=0,
        n_info=0,
    )


def _facts(ctx: AuditContext) -> pd.DataFrame:
    s, g = ctx.balance, ctx.gaps.by_unit
    values = {
        "n_units": s.n_units,
        "n_periods": s.n_periods,
        "n_obs": s.n_obs,
        "T_min": s.t_min,
        "T_mean": s.t_mean,
        "T_max": s.t_max,
        "gamma": s.gamma,
        "nu": s.nu,
        "fill_rate": s.fill_rate,
        "units_with_gaps": int(g["internal_gap"].sum()),
        "units_late_entry": int(g["late_entry"].sum()),
        "units_early_exit": int(g["early_exit"].sum()),
        "units_single_obs": int((ctx.panel.index.counts == 1).sum()),
    }
    if len(ctx.variables) >= 2:
        values["complete_row_share"] = ctx.complete_share
    # a Series keeps the order of the keys on every pandas version (a dict of dicts is sorted by
    # pandas before 2.2)
    return pd.Series(values, name="value").to_frame().rename_axis("fact")


def _variables(ctx: AuditContext) -> pd.DataFrame:
    if not ctx.variables:
        return pd.DataFrame()
    missing = ctx.missing
    assert missing is not None
    table = pd.DataFrame(
        {
            "n_obs": missing.by_variable["n_obs"] - missing.by_variable["n_missing"],
            "pct_missing": missing.by_variable["pct_missing"],
        }
    ).reindex(ctx.variables)
    table["kind"] = ctx.variation["kind"].reindex(ctx.variables) if len(ctx.variation) else None
    table["within_variance_share"] = (
        ctx.within["variance_share"].reindex(ctx.variables) if len(ctx.within) else np.nan
    )
    table["outliers_flagged"] = [
        ctx.outliers(v).n_flagged if v in ctx.observed_variables else 0 for v in ctx.variables
    ]
    table.index.name = "variable"
    return table


@tool(
    "audit",
    ToolSpec(
        summary="The whole panel examined in one call: structure, gaps, missing values, "
        "variation, outliers and equations, with alerts and advice in French and English.",
        validation="alerts raised on panels built for each, tools it aggregates validated "
        "separately, export formats checked",
        similar="expdpy explore_panel_structure; pyxtdpthresh diagnose()",
    ),
)
def audit(
    panel: PanelBase,
    variables: Sequence[str] | None = None,
    *,
    dependent: str | None = None,
    lags: int = 1,
    thresholds: Mapping[str, float] | AuditThresholds | None = None,
    checks: Sequence[str] | None = None,
) -> AuditReport:
    """Examine a declared panel and report what to look at before estimating.

    Parameters
    ----------
    panel
        A declared panel.
    variables
        Variables to examine; by default every numeric variable other than the keys.
    dependent
        Dependent variable of a dynamic model; when given, the audit counts the equations left
        by first differences and by orthogonal deviations with ``lags`` lags.
    lags
        Lags of ``dependent`` in the model.
    thresholds
        Overrides of :class:`~pyxtpanel.explore.checks._context.AuditThresholds`, as a mapping
        (``{"min_units": 50}``) or an instance.
    checks
        Names of the checks to run (default: all of them, see ``CHECKS.names()``).
    """
    names = list(variables) if variables is not None else numeric_names_or_empty(panel)
    limits = (
        thresholds
        if isinstance(thresholds, AuditThresholds)
        else AuditThresholds.from_mapping(thresholds)
    )
    chosen = list(checks) if checks is not None else CHECKS.names()
    for name in chosen:
        if name not in CHECKS:
            raise SpecificationError(
                _i18n.message(
                    "explore.audit.check", name=name, available=", ".join(CHECKS.names())
                ),
                stata_code=198,
                key="explore.audit.check",
            )
    context = AuditContext(panel, names, limits, dependent, lags)
    ordered = sorted(chosen, key=lambda n: (CHECKS.spec(n).order, n))
    found = [alert for name in ordered for alert in run(name, context)]
    rank = {level: i for i, level in enumerate(LEVELS)}
    found.sort(key=lambda a: rank[a.level])  # stable: report order inside a level
    return AuditReport(
        alerts=tuple(found),
        facts=_facts(context),
        variables=_variables(context),
        n_critical=sum(a.level == "critical" for a in found),
        n_warning=sum(a.level == "warning" for a in found),
        n_info=sum(a.level == "info" for a in found),
    )


def numeric_names_or_empty(panel: PanelBase) -> list[str]:
    """Numeric variables other than the keys, or none when there are none."""
    try:
        return numeric_names(panel)
    except SpecificationError:
        return []
