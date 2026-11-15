"""Audit check: values far from the median of their own unit."""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.explore.checks._context import AuditContext
from pyxtpanel.explore.checks._spec import Alert, CheckSpec, Level, check, short_list

_i18n.register(
    "explore",
    {
        "explore.audit.outliers": {
            "en": "{n} values ({share:.2f} % of those scored) of {names} are far from the "
            "median of their unit (|z| > {threshold}).",
            "fr": "{n} valeurs ({share:.2f} % des valeurs évaluées) de {names} sont loin de la "
            "médiane de leur unité (|z| > {threshold}).",
        },
        "explore.audit.outliers.advice": {
            "en": "They are only marked, never removed: inspect them with within_outliers(); "
            "data entry errors are common.",
            "fr": "Elles sont seulement marquées, jamais supprimées : examinez-les avec "
            "within_outliers() ; les erreurs de saisie sont fréquentes.",
        },
    },
)


@check(
    "outliers",
    CheckSpec(
        summary="Within outliers by robust scores.",
        codes=("outliers",),
        order=70,
    ),
)
def outliers(ctx: AuditContext) -> list[Alert]:
    """Raise one alert listing the variables with values far from their unit's median."""
    t = ctx.thresholds
    flagged, scored, names = 0, 0, []
    for name in ctx.variables:
        if name not in ctx.variation.index or ctx.variation.loc[name, "kind"] == "constant":
            continue
        result = ctx.outliers(name)
        if result.n_flagged:
            flagged += result.n_flagged
            scored += result.n_scored
            names.append(f"{name} ({result.n_flagged})")
    if not flagged:
        return []
    share = 100.0 * flagged / max(scored, 1)
    level: Level = "warning" if share >= 100.0 * t.outlier_share else "info"
    return [
        Alert(
            "outliers",
            "outliers",
            level,
            {
                "n": flagged,
                "share": share,
                "names": short_list(names),
                "threshold": t.outlier_threshold,
            },
        )
    ]
