"""Audit check: missing values of the variables."""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.explore.checks._context import AuditContext
from pyxtpanel.explore.checks._spec import Alert, CheckSpec, check, short_list

_i18n.register(
    "explore",
    {
        "explore.audit.all_missing": {
            "en": "No value is observed for: {names}.",
            "fr": "Aucune valeur n'est observée pour : {names}.",
        },
        "explore.audit.all_missing.advice": {
            "en": "Remove these variables or check how they were read.",
            "fr": "Retirez ces variables ou vérifiez comment elles ont été lues.",
        },
        "explore.audit.mostly_missing": {
            "en": "More than {limit:.0f} % missing: {names}.",
            "fr": "Plus de {limit:.0f} % de valeurs manquantes : {names}.",
        },
        "explore.audit.mostly_missing.advice": {
            "en": "A regression using them keeps few rows; check whether they can be rebuilt or "
            "left out.",
            "fr": "Une régression qui les utilise garde peu de lignes ; vérifiez si elles "
            "peuvent être reconstruites ou écartées.",
        },
        "explore.audit.high_missing": {
            "en": "{limit:.0f} % or more missing: {names}.",
            "fr": "{limit:.0f} % ou plus de valeurs manquantes : {names}.",
        },
        "explore.audit.high_missing.advice": {
            "en": "See missing_report() for the periods and units concerned.",
            "fr": "Voir missing_report() pour les périodes et les unités concernées.",
        },
        "explore.audit.listwise_loss": {
            "en": "Keeping the complete rows of the {k} variables leaves {share:.0f} % of the "
            "rows.",
            "fr": "Ne garder que les lignes complètes des {k} variables laisse {share:.0f} % des "
            "lignes.",
        },
        "explore.audit.listwise_loss.advice": {
            "en": "Missing values add up across variables: look at the joint patterns in "
            "missing_report() and at complete_cases().",
            "fr": "Les valeurs manquantes s'additionnent d'une variable à l'autre : regardez les "
            "motifs joints de missing_report() et complete_cases().",
        },
    },
)


@check(
    "missing",
    CheckSpec(
        summary="Missing values by variable and in the complete rows.",
        codes=("all_missing", "mostly_missing", "high_missing", "listwise_loss"),
        order=40,
    ),
)
def missing(ctx: AuditContext) -> list[Alert]:
    """Raise alerts about missing values."""
    report = ctx.missing
    if report is None:
        return []
    t = ctx.thresholds
    share = report.by_variable["pct_missing"] / 100.0
    everything = [str(v) for v in share.index[share >= 1.0]]
    heavy = [str(v) for v in share.index[(share >= t.missing_critical) & (share < 1.0)]]
    some = [str(v) for v in share.index[(share >= t.missing_share) & (share < t.missing_critical)]]
    alerts = []
    if everything:
        alerts.append(
            Alert("missing", "all_missing", "critical", {"names": short_list(everything)})
        )
    if heavy:
        alerts.append(
            Alert(
                "missing",
                "mostly_missing",
                "critical",
                {"names": short_list(heavy), "limit": 100.0 * t.missing_critical},
            )
        )
    if some:
        alerts.append(
            Alert(
                "missing",
                "high_missing",
                "warning",
                {"names": short_list(some), "limit": 100.0 * t.missing_share},
            )
        )
    if len(ctx.variables) >= 2 and ctx.complete_share < t.complete_share:
        alerts.append(
            Alert(
                "missing",
                "listwise_loss",
                "warning",
                {"k": len(ctx.variables), "share": 100.0 * ctx.complete_share},
            )
        )
    return alerts
