"""Audit check: where the variation of each variable comes from."""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.explore.checks._context import AuditContext
from pyxtpanel.explore.checks._spec import Alert, CheckSpec, Level, check, short_list

_i18n.register(
    "explore",
    {
        "explore.audit.constant_variable": {
            "en": "Constant in the whole panel: {names}.",
            "fr": "Constante dans tout le panel : {names}.",
        },
        "explore.audit.constant_variable.advice": {
            "en": "A constant carries no information; drop it (the model has its own intercept).",
            "fr": "Une constante n'apporte rien ; retirez-la (le modèle a sa propre constante).",
        },
        "explore.audit.time_invariant": {
            "en": "Varies only between units: {names}.",
            "fr": "Ne varie qu'entre les unités : {names}.",
        },
        "explore.audit.time_invariant.advice": {
            "en": "Absorbed by unit fixed effects: its coefficient is not identified in a "
            "within or first-difference model.",
            "fr": "Absorbée par les effets fixes d'unité : son coefficient n'est pas identifié "
            "dans un modèle within ou en différences premières.",
        },
        "explore.audit.common_variable": {
            "en": "Identical for all units in each period: {names}.",
            "fr": "Identique pour toutes les unités à chaque période : {names}.",
        },
        "explore.audit.common_variable.advice": {
            "en": "Absorbed by period effects; keep it only if the model has no time dummies.",
            "fr": "Absorbée par les effets de période ; ne la gardez que si le modèle n'a pas "
            "d'indicatrices de temps.",
        },
        "explore.audit.little_within": {
            "en": "Less than {limit:.0f} % of the variance is left after removing the unit "
            "means: {names}.",
            "fr": "Moins de {limit:.0f} % de la variance reste une fois les moyennes d'unité "
            "retirées : {names}.",
        },
        "explore.audit.little_within.advice": {
            "en": "Within and first-difference estimators will be imprecise for these "
            "variables; see residual_variation().",
            "fr": "Les estimateurs within et en différences premières seront imprécis pour ces "
            "variables ; voir residual_variation().",
        },
    },
)


@check(
    "variation",
    CheckSpec(
        summary="Constant, time-invariant, common and almost time-invariant variables.",
        codes=("constant_variable", "time_invariant", "common_variable", "little_within"),
        order=60,
    ),
)
def variation(ctx: AuditContext) -> list[Alert]:
    """Raise alerts about the variation of the variables."""
    table = ctx.variation
    if table.empty:
        return []
    kind = table["kind"]
    alerts = []
    rules: tuple[tuple[str, str, Level], ...] = (
        ("constant", "constant_variable", "critical"),
        ("time_invariant", "time_invariant", "warning"),
        ("common", "common_variable", "info"),
    )
    for label, code, level in rules:
        names = [str(v) for v in table.index[kind == label]]
        if names:
            alerts.append(Alert("variation", code, level, {"names": short_list(names)}))
    share = ctx.within["variance_share"]
    limit = ctx.thresholds.within_share
    low = [str(v) for v in table.index[kind == "varying"] if v in share.index and share[v] < limit]
    if low:
        alerts.append(
            Alert(
                "variation",
                "little_within",
                "warning",
                {"names": short_list(low), "limit": 100.0 * limit},
            )
        )
    return alerts
