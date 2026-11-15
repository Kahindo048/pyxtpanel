"""Statistical tests of panel data, categories A (on data) and B (generic post-estimation).

Layer 4.  Every test is a pure function that returns a :class:`TestResult`; the estimation
packages call the same functions as users do, with bare arrays for the category B tests
(``from pyxtpanel.diagnostics import hansen, ar_test``).  This import path is the stable surface
the estimators depend on.
"""

from pyxtpanel.diagnostics._core import (
    Equations,
    Note,
    TestResult,
    TestSpec,
    list_tests,
    panel_tests,
    register_test,
    resolve_test,
    run_panel_test,
    run_test,
)
from pyxtpanel.diagnostics._shared.overid import moment_covariance
from pyxtpanel.diagnostics.crosssection.cross_correlation import CrossCorrelation, cross_correlation
from pyxtpanel.diagnostics.crosssection.frees import frees
from pyxtpanel.diagnostics.crosssection.friedman import friedman
from pyxtpanel.diagnostics.crosssection.lm_cross_section import lm_cross_section
from pyxtpanel.diagnostics.crosssection.pesaran_cd import pesaran_cd
from pyxtpanel.diagnostics.effects.f_effects import f_effects
from pyxtpanel.diagnostics.effects.lm_random_effects import lm_random_effects
from pyxtpanel.diagnostics.heteroskedasticity.breusch_pagan import breusch_pagan
from pyxtpanel.diagnostics.heteroskedasticity.wald_groupwise import wald_groupwise
from pyxtpanel.diagnostics.heteroskedasticity.white import white
from pyxtpanel.diagnostics.normality.jarque_bera import jarque_bera
from pyxtpanel.diagnostics.normality.shapiro_francia import shapiro_francia
from pyxtpanel.diagnostics.normality.shapiro_wilk import shapiro_wilk
from pyxtpanel.diagnostics.normality.skewness_kurtosis import skewness_kurtosis
from pyxtpanel.diagnostics.postestimation.ar_test import ar_test
from pyxtpanel.diagnostics.postestimation.diff_hansen import diff_hansen
from pyxtpanel.diagnostics.postestimation.hansen import hansen
from pyxtpanel.diagnostics.postestimation.hausman import hausman
from pyxtpanel.diagnostics.postestimation.sargan import sargan
from pyxtpanel.diagnostics.postestimation.wald import wald
from pyxtpanel.diagnostics.serial.wooldridge import wooldridge
from pyxtpanel.diagnostics.specification.reset import reset
from pyxtpanel.diagnostics.specification.vif import VarianceInflation, vif
from pyxtpanel.diagnostics.unitroot.fisher import fisher
from pyxtpanel.diagnostics.unitroot.harris_tzavalis import harris_tzavalis
from pyxtpanel.diagnostics.unitroot.im_pesaran_shin import im_pesaran_shin
from pyxtpanel.diagnostics.unitroot.levin_lin_chu import levin_lin_chu

__all__ = [
    "CrossCorrelation",
    "VarianceInflation",
    "cross_correlation",
    "f_effects",
    "frees",
    "friedman",
    "fisher",
    "harris_tzavalis",
    "im_pesaran_shin",
    "jarque_bera",
    "levin_lin_chu",
    "lm_cross_section",
    "lm_random_effects",
    "pesaran_cd",
    "shapiro_francia",
    "shapiro_wilk",
    "skewness_kurtosis",
    "wald_groupwise",
    "Equations",
    "Note",
    "TestResult",
    "TestSpec",
    "ar_test",
    "breusch_pagan",
    "diff_hansen",
    "hansen",
    "hausman",
    "list_tests",
    "moment_covariance",
    "panel_tests",
    "register_test",
    "reset",
    "resolve_test",
    "run_panel_test",
    "run_test",
    "sargan",
    "vif",
    "wald",
    "white",
    "wooldridge",
]
