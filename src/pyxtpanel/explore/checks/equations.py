"""Audit check: equations left for a dynamic estimator (needs a dependent variable)."""

from __future__ import annotations

from pyxtpanel import _i18n
from pyxtpanel.explore.checks._context import AuditContext
from pyxtpanel.explore.checks._spec import Alert, CheckSpec, check

_i18n.register(
    "explore",
    {
        "explore.audit.no_equations": {
            "en": "No usable equation for {variable} with {lags} lag(s), in first differences "
            "or in orthogonal deviations.",
            "fr": "Aucune équation utilisable pour {variable} avec {lags} retard(s), ni en "
            "différences premières ni en écarts orthogonaux.",
        },
        "explore.audit.no_equations.advice": {
            "en": "Reduce the number of lags or regressors: the units are too short or have "
            "too many gaps. See usable_equations().",
            "fr": "Réduisez le nombre de retards ou de régresseurs : les unités sont trop "
            "courtes ou ont trop de trous. Voir usable_equations().",
        },
        "explore.audit.few_equations": {
            "en": "Only {share:.0f} % of the observations of {variable} give a first-difference "
            "equation ({usable} of {observed}).",
            "fr": "Seulement {share:.0f} % des observations de {variable} donnent une équation "
            "en différences premières ({usable} sur {observed}).",
        },
        "explore.audit.few_equations.advice": {
            "en": "Gaps and short units remove equations; see usable_equations().",
            "fr": "Les trous et les unités courtes retirent des équations ; voir "
            "usable_equations().",
        },
        "explore.audit.prefer_fod": {
            "en": "Orthogonal deviations keep {fod} equations against {fd} for first differences.",
            "fr": "Les écarts orthogonaux conservent {fod} équations contre {fd} en différences "
            "premières.",
        },
        "explore.audit.prefer_fod.advice": {
            "en": "Prefer forward orthogonal deviations: this panel has gaps.",
            "fr": "Préférez les écarts orthogonaux vers l'avant : ce panel a des trous.",
        },
    },
)


@check(
    "equations",
    CheckSpec(
        summary="Equations left by first differences and orthogonal deviations.",
        codes=("no_equations", "few_equations", "prefer_fod"),
        order=80,
    ),
)
def equations(ctx: AuditContext) -> list[Alert]:
    """Raise alerts about the equations a dynamic estimator can use (when ``dependent`` is set)."""
    if ctx.dependent is None:
        return []
    fd, fod = ctx.equations("fd"), ctx.equations("fod")
    base = {"variable": ctx.dependent, "lags": ctx.lags}
    if fd.n_usable == 0 and fod.n_usable == 0:
        return [Alert("equations", "no_equations", "critical", base)]
    alerts = []
    observed = max(fd.n_observed, 1)
    if fd.n_usable / observed < ctx.thresholds.equations_share:
        alerts.append(
            Alert(
                "equations",
                "few_equations",
                "warning",
                {
                    **base,
                    "share": 100.0 * fd.n_usable / observed,
                    "usable": fd.n_usable,
                    "observed": fd.n_observed,
                },
            )
        )
    if fod.n_usable > fd.n_usable * (1.0 + ctx.thresholds.fod_gain):
        alerts.append(
            Alert(
                "equations",
                "prefer_fod",
                "info",
                {"fod": fod.n_usable, "fd": fd.n_usable},
            )
        )
    return alerts
