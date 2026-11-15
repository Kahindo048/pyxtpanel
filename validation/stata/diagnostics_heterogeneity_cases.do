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
mata: pyx_export("../out/sh_csa__sh.csv", (st_numscalar("d1") \ st_numscalar("d2") \ st_numscalar("p1") \ st_numscalar("p2") \ st_numscalar("sw") \ st_numscalar("swp")))
set type float

exit, clear
