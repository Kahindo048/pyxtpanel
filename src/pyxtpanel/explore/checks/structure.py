"""Audit check: size of the panel (units, periods, units with one observation)."""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.explore.checks._context import AuditContext
from pyxtpanel.explore.checks._spec import Alert, CheckSpec, Level, check

_i18n.register(
    "explore",
    {
        "explore.audit.single_unit": {
            "en": "The panel has a single unit: there is no cross-section.",
            "fr": "Le panel n'a qu'une seule unité : il n'y a pas de coupe transversale.",
        },
        "explore.audit.single_unit.advice": {
            "en": "Check the identifier; panel methods need several units.",
            "fr": "Vérifiez l'identifiant ; les méthodes de panel exigent plusieurs unités.",
        },
        "explore.audit.few_units": {
            "en": "Only {n} units (fewer than {limit}): standard errors and between-unit "
            "estimates are fragile.",
            "fr": "Seulement {n} unités (moins de {limit}) : les écarts types et les estimations "
            "between sont fragiles.",
        },
        "explore.audit.few_units.advice": {
            "en": "Use small-sample corrections or a bootstrap, and be careful with clustering "
            "by unit.",
            "fr": "Utilisez des corrections de petit échantillon ou un bootstrap, et soyez "
            "prudent avec le regroupement (cluster) par unité.",
        },
        "explore.audit.very_short": {
            "en": "The calendar has {t} period(s): lags and first differences cannot be formed.",
            "fr": "Le calendrier compte {t} période(s) : retards et différences premières "
            "impossibles.",
        },
        "explore.audit.very_short.advice": {
            "en": "A dynamic or first-difference model needs at least 3 periods.",
            "fr": "Un modèle dynamique ou en différences premières exige au moins 3 périodes.",
        },
        "explore.audit.short_panel": {
            "en": "Short panel: {t} periods (fewer than {limit}).",
            "fr": "Panel court : {t} périodes (moins de {limit}).",
        },
        "explore.audit.short_panel.advice": {
            "en": "Dynamic estimators with many instrument lags are weak here: limit the lags "
            "and collapse the instruments.",
            "fr": "Les estimateurs dynamiques avec beaucoup de retards d'instruments sont "
            "faibles ici : limitez les retards et réduisez les instruments (collapse).",
        },
        "explore.audit.single_obs_units": {
            "en": "{n} units ({share:.0f} %) have a single observation: they carry no within "
            "information.",
            "fr": "{n} unités ({share:.0f} %) n'ont qu'une observation : elles n'apportent "
            "aucune information within.",
        },
        "explore.audit.single_obs_units.advice": {
            "en": "They only contribute to between estimates; consider leaving them out of "
            "within and dynamic models.",
            "fr": "Elles ne servent qu'aux estimations between ; envisagez de les écarter des "
            "modèles within et dynamiques.",
        },
    },
)


@check(
    "structure",
    CheckSpec(
        summary="Number of units and periods, units observed once.",
        codes=("single_unit", "few_units", "very_short", "short_panel", "single_obs_units"),
        order=10,
    ),
)
def structure(ctx: AuditContext) -> list[Alert]:
    """Raise alerts about the size of the panel."""
    t, s = ctx.thresholds, ctx.balance
    alerts = []
    if s.n_units == 1:
        alerts.append(Alert("structure", "single_unit", "critical"))
    elif s.n_units < t.min_units:
        alerts.append(
            Alert("structure", "few_units", "warning", {"n": s.n_units, "limit": t.min_units})
        )
    if s.n_periods < 3:
        alerts.append(Alert("structure", "very_short", "critical", {"t": s.n_periods}))
    elif s.n_periods < t.min_periods:
        alerts.append(
            Alert("structure", "short_panel", "warning", {"t": s.n_periods, "limit": t.min_periods})
        )
    single = int((ctx.panel.index.counts == 1).sum())
    if single and s.n_units > 1:
        share = 100.0 * single / s.n_units
        level: Level = "warning" if single / s.n_units > t.single_obs_share else "info"
        alerts.append(Alert("structure", "single_obs_units", level, {"n": single, "share": share}))
    return alerts
