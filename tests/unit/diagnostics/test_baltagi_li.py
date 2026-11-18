"""Baltagi-Li tests: R's plm values, an oracle in plain loops, hand-computed cases, refusals."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from pyxtpanel._errors import DataError, NumericalError, SpecificationError
from pyxtpanel.diagnostics import baltagi_li, run_test

KEYS = (
    "joint",
    "re",
    "re_adjusted",
    "serial",
    "serial_adjusted",
    "re_onesided",
    "re_adjusted_onesided",
)


def balanced_residuals() -> tuple[np.ndarray, np.ndarray]:
    """Pooled residuals of the panel (30 x 8, unit effect and AR(0.4) errors) given to plm."""
    rng = np.random.default_rng(11)
    n, t = 30, 8
    u = np.repeat(rng.normal(size=n) * 0.8, t)
    v = np.zeros((n, t))
    for period in range(t):
        v[:, period] = (0.4 * v[:, period - 1] if period else 0) + rng.normal(size=n)
    x = rng.normal(size=n * t)
    y = 1 + 0.5 * x + u + v.ravel()
    design = np.column_stack([np.ones(n * t), x])
    e = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return e, np.arange(0, n * t + 1, t)


def oracle(e: np.ndarray, offsets: np.ndarray) -> dict[str, float]:
    """The seven statistics from the definitions, unit by unit in plain loops."""
    m, groups = len(e), len(offsets) - 1
    sum_sq = sum(v * v for v in e)
    sums = a = num = den = 0.0
    for g in range(groups):
        unit = [float(v) for v in e[offsets[g] : offsets[g + 1]]]
        sums += sum(unit) ** 2
        a += len(unit) ** 2
        for k in range(1, len(unit)):
            num += unit[k] * unit[k - 1]
            den += unit[k] ** 2
    a1, b1 = 1 - sums / sum_sq, num / den
    lm_re = m**2 * a1**2 / (2 * (a - m))
    lm_se = m**2 * b1**2 / (m - groups)
    alm_re = m**2 * (a1 + 2 * b1) ** 2 / (2 * (a - 3 * m + 2 * groups))
    alm_se = lm_se - lm_re + alm_re
    return {
        "joint": lm_re + alm_se,
        "re": lm_re,
        "re_adjusted": alm_re,
        "serial": lm_se,
        "serial_adjusted": alm_se,
        "re_onesided": -math.sqrt(m**2 / (2 * (a - m))) * a1,
        "re_adjusted_onesided": -math.sqrt(m**2 / (2 * (a - 3 * m + 2 * groups))) * (a1 + 2 * b1),
    }


# plm 2.6.7, R 4.6.1 on the panel above: pbsytest(y ~ x, test = "j" | "re" | "ar") and
# plmtest(y ~ x, type = "bp"); plm's "re" is the one-sided standard normal statistic and its "ar"
# the chi-square(1) robust statistic of serial correlation.
PLM = {
    "joint": (194.98175944752, 4.57350228158176e-43),
    "re_adjusted_onesided": (10.9278356960991, 4.24239851431702e-28),
    "serial_adjusted": (4.26107697480489, 0.0389951917085724),
    "re": (190.720682472715, None),
}


@pytest.mark.parametrize("key", list(PLM))
def test_the_statistics_are_those_of_plm(key: str) -> None:
    e, offsets = balanced_residuals()
    result = baltagi_li(e, offsets, statistic=key)
    statistic, pvalue = PLM[key]
    assert result.statistic == pytest.approx(statistic, rel=1e-12)
    if pvalue is not None:
        assert result.pvalue == pytest.approx(pvalue, rel=1e-10)


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_an_unbalanced_panel_agrees_with_the_oracle_in_loops(seed: int) -> None:
    rng = np.random.default_rng(seed)
    sizes = rng.integers(3, 9, size=25)
    offsets = np.concatenate([[0], np.cumsum(sizes)])
    e = rng.normal(size=int(offsets[-1])) + np.repeat(rng.normal(size=25), sizes)
    want = oracle(e, offsets)
    for key in KEYS:
        assert baltagi_li(e, offsets, statistic=key).statistic == pytest.approx(
            want[key], rel=1e-12
        )


def test_the_joint_statistic_is_the_marginal_one_plus_the_robust_serial_one() -> None:
    e, offsets = balanced_residuals()
    value = {k: baltagi_li(e, offsets, statistic=k).statistic for k in KEYS}
    assert value["joint"] == pytest.approx(value["re"] + value["serial_adjusted"], rel=1e-13)
    assert value["serial_adjusted"] == pytest.approx(
        value["serial"] - value["re"] + value["re_adjusted"], rel=1e-13
    )


def test_the_degrees_of_freedom_and_tails_follow_the_distribution_of_each_statistic() -> None:
    e, offsets = balanced_residuals()
    joint = baltagi_li(e, offsets)
    assert joint.df == (2.0,) and joint.statistic_name == "chi2"
    assert joint.pvalue == pytest.approx(float(stats.chi2.sf(joint.statistic, 2)), rel=1e-13)
    one_sided = baltagi_li(e, offsets, statistic="re_onesided")
    assert one_sided.df == () and one_sided.statistic_name == "z"
    assert one_sided.pvalue == pytest.approx(float(stats.norm.sf(one_sided.statistic)), rel=1e-13)
    assert list(joint.details.index) == list(KEYS) and joint.n_groups == 30 and joint.n_obs == 240


def test_hand_computed_units_of_two_observations() -> None:
    # e = (1, 1 | -1, -1): A = 1 - 8/4 = -1, B = (1 + 1) / (1 + 1) = 1, m = 4, G = 2, S = 8
    e, offsets = [1.0, 1.0, -1.0, -1.0], [0, 2, 4]
    assert baltagi_li(e, offsets, statistic="re").statistic == pytest.approx(2.0, rel=1e-14)
    assert baltagi_li(e, offsets, statistic="serial").statistic == pytest.approx(8.0, rel=1e-14)
    # no unit has three observations: the robust statistics divide by zero and are refused
    for key in ("joint", "re_adjusted", "serial_adjusted", "re_adjusted_onesided"):
        with pytest.raises(NumericalError, match="three observations"):
            baltagi_li(e, offsets, statistic=key)


def test_units_of_one_observation_leave_no_serial_statistic() -> None:
    with pytest.raises(NumericalError, match="two observations"):
        baltagi_li([1.0, 2.0, 3.0], [0, 1, 2, 3], statistic="serial")


def test_it_tells_random_effects_from_serial_correlation() -> None:
    rng = np.random.default_rng(5)
    offsets = np.arange(0, 401, 8)
    effect = np.repeat(rng.normal(size=50), 8)
    # random effects alone: the robust serial test does not reject, the robust effect test does
    e = effect + rng.normal(size=400)
    assert baltagi_li(e, offsets, statistic="re_adjusted").rejects(0.001)
    assert not baltagi_li(e, offsets, statistic="serial_adjusted").rejects(0.001)
    # AR(1) errors alone: the reverse
    v = np.zeros((50, 8))
    for t in range(8):
        v[:, t] = (0.6 * v[:, t - 1] if t else 0) + rng.normal(size=50)
    assert baltagi_li(v.ravel(), offsets, statistic="serial_adjusted").rejects(0.001)


def test_a_unit_that_skips_a_period_is_refused_with_the_code_of_xttest1() -> None:
    e = np.arange(1.0, 8.0) ** 0.5
    time = np.array([1, 2, 4, 5, 1, 2, 3])
    with pytest.raises(DataError, match="has a gap") as refusal:
        baltagi_li(e, [0, 4, 7], time_index=time)
    assert refusal.value.stata_code == 498
    assert baltagi_li(e, [0, 4, 7], time_index=[1, 2, 3, 4, 1, 2, 3]).n_groups == 2


def test_refusals() -> None:
    e = np.arange(1.0, 7.0)
    with pytest.raises(SpecificationError, match="statistic must be one of"):
        baltagi_li(e, [0, 3, 6], statistic="lm")
    with pytest.raises(SpecificationError, match="offsets must start at 0"):
        baltagi_li(e, [0, 3])
    with pytest.raises(SpecificationError, match="one entry per residual"):
        baltagi_li(e, [0, 3, 6], time_index=[1, 2, 3])
    with pytest.raises(NumericalError, match="all residuals are zero"):
        baltagi_li(np.zeros(6), [0, 3, 6])


def test_by_name_and_alias_in_both_languages() -> None:
    e, offsets = balanced_residuals()
    result = run_test("xttest1", e, offsets)
    assert result.name == "baltagi_li"
    assert "Baltagi-Li tests" in result.summary("en")
    assert "Tests de Baltagi-Li" in result.summary("fr")
