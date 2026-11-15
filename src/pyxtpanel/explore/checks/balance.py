"""Audit check: balance, gaps, entries and exits of the units."""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.explore.checks._context import AuditContext
from pyxtpanel.explore.checks._spec import Alert, CheckSpec, Level, check

_i18n.register(
    "explore",
    {
        "explore.audit.unbalanced": {
            "en": "The panel is unbalanced: gamma = {gamma:.2f}, nu = {nu:.2f}, {fill:.0f} % of "
            "the unit-period cells observed.",
            "fr": "Le panel est déséquilibré : gamma = {gamma:.2f}, nu = {nu:.2f}, {fill:.0f} % "
            "des cellules unité-période observées.",
        },
        "explore.audit.unbalanced.advice": {
            "en": "Compare the ways of balancing with balance_cost() before choosing; most "
            "estimators accept an unbalanced panel.",
            "fr": "Comparez les façons de cylindrer avec balance_cost() avant de choisir ; la "
            "plupart des estimateurs acceptent un panel déséquilibré.",
        },
        "explore.audit.internal_gaps": {
            "en": "{n} units ({share:.0f} %) have gaps inside their span ({periods} periods "
            "missing in all).",
            "fr": "{n} unités ({share:.0f} %) ont des trous dans leur durée ({periods} périodes "
            "manquantes au total).",
        },
        "explore.audit.internal_gaps.advice": {
            "en": "First differences lose equations at every gap; forward orthogonal deviations "
            "keep them. See usable_equations() and gap_types().",
            "fr": "Les différences premières perdent des équations à chaque trou ; les écarts "
            "orthogonaux les conservent. Voir usable_equations() et gap_types().",
        },
        "explore.audit.heavy_attrition": {
            "en": "{n} units ({share:.0f} %) stop being observed before the last period.",
            "fr": "{n} unités ({share:.0f} %) cessent d'être observées avant la dernière période.",
        },
        "explore.audit.heavy_attrition.advice": {
            "en": "Check whether leaving is related to the outcome (selection); compare the "
            "units that stay with those that leave (attrition()).",
            "fr": "Vérifiez si la sortie est liée à la variable d'intérêt (sélection) ; "
            "comparez les unités qui restent à celles qui partent (attrition()).",
        },
        "explore.audit.late_entry": {
            "en": "{n} units ({share:.0f} %) are first observed after the first period.",
            "fr": "{n} unités ({share:.0f} %) sont observées pour la première fois après la "
            "première période.",
        },
        "explore.audit.late_entry.advice": {
            "en": "Entry may be informative (new firms, new members); use an estimator that "
            "accepts an unbalanced panel.",
            "fr": "L'entrée peut être informative (nouvelles entreprises, nouveaux membres) ; "
            "utilisez un estimateur qui accepte un panel déséquilibré.",
        },
    },
)


@check(
    "balance",
    CheckSpec(
        summary="Imbalance measures, internal gaps, late entries and early exits.",
        codes=("unbalanced", "internal_gaps", "heavy_attrition", "late_entry"),
        order=20,
    ),
)
def balance(ctx: AuditContext) -> list[Alert]:
    """Raise alerts about balance, gaps and attrition."""
    t, s, g = ctx.thresholds, ctx.balance, ctx.gaps
    n = s.n_units
    alerts = []
    if s.gamma < t.gamma or s.nu < t.gamma:
        alerts.append(
            Alert(
                "balance",
                "unbalanced",
                "warning",
                {"gamma": s.gamma, "nu": s.nu, "fill": 100.0 * s.fill_rate},
            )
        )
    with_gap = int(g.by_unit["internal_gap"].sum())
    if with_gap:
        level: Level = "warning" if with_gap / n >= t.gap_share else "info"
        alerts.append(
            Alert(
                "balance",
                "internal_gaps",
                level,
                {
                    "n": with_gap,
                    "share": 100.0 * with_gap / n,
                    "periods": int(g.by_unit["gap_periods"].sum()),
                },
            )
        )
    leaving = int(g.by_unit["early_exit"].sum())
    if leaving and leaving / n >= t.exit_share:
        alerts.append(
            Alert(
                "balance",
                "heavy_attrition",
                "warning",
                {"n": leaving, "share": 100.0 * leaving / n},
            )
        )
    entering = int(g.by_unit["late_entry"].sum())
    if entering and entering / n >= t.exit_share:
        alerts.append(
            Alert(
                "balance",
                "late_entry",
                "info",
                {"n": entering, "share": 100.0 * entering / n},
            )
        )
    return alerts
