version 17
set more off
clear all

mata:
void pyx_export(string scalar fname, real matrix M)
{
    real scalar fh, i, j
    string scalar line
    if (fileexists(fname)) unlink(fname)
    fh = fopen(fname, "w")
    for (i = 1; i <= rows(M); i++) {
        line = strofreal(M[i,1], "%21.15g")
        for (j = 2; j <= cols(M); j++) line = line + "," + strofreal(M[i,j], "%21.15g")
        fput(fh, line)
    }
    fclose(fh)
}
end

mata: pyx_export("../out/_meta__stata.csv", (st_numscalar("c(stata_version)") \ st_numscalar("c(MP)")))

* ---- case sh_homogeneous
use "../data/sim_slopes_homog.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_homogeneous__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_heterogeneous
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_heterogeneous__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_unbalanced
use "../data/sim_slopes_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_unbalanced__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_partial
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, partial(w) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_partial__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_noconstant
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, noconstant nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_noconstant__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_ar
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, ar nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_ar__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_csa
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, crosssectional(y x z) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_csa__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_csa_lags
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, crosssectional(y x z, cr_lags(2)) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_csa_lags__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

* ---- case sh_hac_bartlett
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, hac nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_bartlett__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_qs
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, hac kernel(qs) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_qs__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_truncated
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, hac kernel(truncated) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_truncated__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_whitening
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, whitening nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_whitening__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_bandwidth
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, bw(3) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_bandwidth__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_qs_bandwidth
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, kernel(qs) bw(4) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_qs_bandwidth__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_truncated_whitening
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, whitening kernel(truncated) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_truncated_whitening__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_partial
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, hac partial(w) crosssectional(y x z) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_partial__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_unbalanced
use "../data/sim_slopes_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, hac nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_unbalanced__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_unbalanced_qs
use "../data/sim_slopes_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, hac kernel(qs) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_unbalanced_qs__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_unbalanced_whitening
use "../data/sim_slopes_unbalanced.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, whitening kernel(truncated) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_unbalanced_whitening__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_bandwidth_zero
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, bw(0) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_bandwidth_zero__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_csa_lags
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z, hac crosssectional(y x z, cr_lags(1)) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (2 + d1 * sqrt(2 * 2) / sqrt(nunits))
scalar swp = chi2tail(2 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_csa_lags__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_one_regressor
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x, hac nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (1 + d1 * sqrt(2 * 1) / sqrt(nunits))
scalar swp = chi2tail(1 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_one_regressor__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_one_whitening
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x, whitening kernel(qs) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (1 + d1 * sqrt(2 * 1) / sqrt(nunits))
scalar swp = chi2tail(1 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_one_whitening__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_one_truncated
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x, hac kernel(truncated) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (1 + d1 * sqrt(2 * 1) / sqrt(nunits))
scalar swp = chi2tail(1 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_one_truncated__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_three_regressors
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z w, hac nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (3 + d1 * sqrt(2 * 3) / sqrt(nunits))
scalar swp = chi2tail(3 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_three_regressors__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

* ---- case sh_hac_three_truncated
use "../data/sim_slopes_het.dta", clear
quietly xtset id t
set type double
quietly xtdescribe
scalar nunits = r(N)
quietly xthst y x z w, hac kernel(truncated) nooutput
matrix dd = r(delta)
scalar d1 = dd[1, 1]
scalar d2 = dd[2, 1]
scalar p1 = 2 * normal(-abs(d1))
scalar p2 = 2 * normal(-abs(d2))
scalar sw = nunits * (3 + d1 * sqrt(2 * 3) / sqrt(nunits))
scalar swp = chi2tail(3 * (nunits - 1), sw)
capture scalar bw = r(bw)
mata: pyx_export("../out/sh_hac_three_truncated__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp") \ st_numscalar("bw")))
set type float

exit, clear
