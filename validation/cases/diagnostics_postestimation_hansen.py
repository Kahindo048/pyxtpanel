"""Hansen J against ``ivregress gmm`` + ``estat overid`` in Stata 17 (domain
``diagnostics_postestimation``)."""

from __future__ import annotations

import numpy as np

from pyxtpanel.diagnostics.postestimation.hansen import hansen
from validation.cases._iv import fit
from validation.registry import case, export_scalars

DOMAIN = "diagnostics_postestimation"


def _python(panel: str, instruments: tuple[str, ...]) -> dict[str, np.ndarray]:
    f = fit(panel, instruments)
    result = hansen(f.hansen_moments, f.hansen_weight, f.n_parameters, n_groups=f.n_groups)
    return {"j": np.array([result.statistic, result.df[0], result.pvalue])}


def _stata(tag: str, instruments: str) -> list[str]:
    return [
        "quietly xtset id t",
        f"quietly ivregress gmm y (x = {instruments}), wmatrix(cluster id)",
        "quietly estat overid",
        "scalar hj = r(HansenJ)",
        "scalar hdf = r(J_df)",
        "scalar hp = r(p_HansenJ)",
        export_scalars(tag, "j", ["hj", "hdf", "hp"]),
    ]


@case(
    "hansen_ivregress",
    domain=DOMAIN,
    dataset="sim_iv_valid",
    stata=_stata("hansen_ivregress", "z1 z2 z3 z4"),
    atol=1e-9,
    rtol=1e-9,
)
def _valid() -> dict[str, np.ndarray]:
    return _python("iv_valid", ("z1", "z2", "z3", "z4"))


@case(
    "hansen_ivregress_invalid",
    domain=DOMAIN,
    dataset="sim_iv_invalid",
    stata=_stata("hansen_ivregress_invalid", "z1 z2 z3 z4"),
    atol=1e-9,
    rtol=1e-9,
)
def _invalid() -> dict[str, np.ndarray]:
    return _python("iv_invalid", ("z1", "z2", "z3", "z4"))


@case(
    "hansen_just_identified",
    domain=DOMAIN,
    dataset="sim_iv_valid",
    stata=[
        "quietly xtset id t",
        "quietly ivregress gmm y (x = z1), wmatrix(cluster id)",
        "capture estat overid",
        "scalar refused = (_rc != 0)",
        export_scalars("hansen_just_identified", "refused", ["refused"]),
    ],
)
def _just_identified() -> dict[str, np.ndarray]:
    f = fit("iv_valid", ("z1",))
    result = hansen(f.hansen_moments, f.hansen_weight, f.n_parameters)
    return {"refused": np.array([float(np.isnan(result.statistic))])}
