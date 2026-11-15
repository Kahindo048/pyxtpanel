"""Sargan against ``ivregress 2sls`` + ``estat overid`` in Stata 17 (domain
``diagnostics_postestimation``)."""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.postestimation.sargan import sargan
from validation.cases._iv import fit
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_postestimation"


def _python(panel: str, instruments: tuple[str, ...]) -> dict[str, np.ndarray]:
    f = fit(panel, instruments)
    result = sargan(
        f.sargan_moments,
        f.sargan_weight,
        f.n_parameters,
        scale=f.sargan_scale,
        n_groups=f.n_groups,
    )
    return {"s": np.array([result.statistic, result.df[0], result.pvalue])}


def _stata(tag: str, instruments: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly ivregress 2sls y (x = {instruments})",
        "quietly estat overid",
        "scalar ss = r(sargan)",
        "scalar sdf = r(df)",
        "scalar sp = r(p_sargan)",
        export_scalars(tag, "s", ["ss", "sdf", "sp"]),
    ]


@case(
    "sargan_ivregress",
    domain=DOMAIN,
    dataset="sim_iv_valid",
    stata=_stata("sargan_ivregress", "z1 z2 z3 z4"),
    atol=1e-9,
    rtol=1e-9,
)
def _valid() -> dict[str, np.ndarray]:
    return _python("iv_valid", ("z1", "z2", "z3", "z4"))


@case(
    "sargan_ivregress_invalid",
    domain=DOMAIN,
    dataset="sim_iv_invalid",
    stata=_stata("sargan_ivregress_invalid", "z1 z2 z3 z4"),
    atol=1e-9,
    rtol=1e-9,
)
def _invalid() -> dict[str, np.ndarray]:
    return _python("iv_invalid", ("z1", "z2", "z3", "z4"))


@case(
    "sargan_just_identified",
    domain=DOMAIN,
    dataset="sim_iv_valid",
    stata=[
        "quietly xtset id t",
        "quietly ivregress 2sls y (x = z1)",
        "capture estat overid",
        "scalar refused = (_rc != 0)",
        export_scalars("sargan_just_identified", "refused", ["refused"]),
    ],
)
def _just_identified() -> dict[str, np.ndarray]:
    f = fit("iv_valid", ("z1",))
    result = sargan(f.sargan_moments, f.sargan_weight, f.n_parameters, scale=f.sargan_scale)
    return {"refused": np.array([float(np.isnan(result.statistic))])}
