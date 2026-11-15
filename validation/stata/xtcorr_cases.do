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

* ---- case xtcorr_balanced
use "../data/sim_balanced.dta", clear
xtset id t
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_balanced__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_balanced__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

* ---- case xtcorr_gaps
use "../data/sim_gaps.dta", clear
xtset id t
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_gaps__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_gaps__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

* ---- case xtcorr_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_entry_exit__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_entry_exit__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

* ---- case xtcorr_weak
use "../data/sim_weak.dta", clear
xtset id t
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_weak__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_weak__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

* ---- case xtcorr_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_delta2__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_delta2__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

* ---- case xtcorr_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_unsorted_ids__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_unsorted_ids__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

* ---- case xtcorr_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_one_obs_unit__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_one_obs_unit__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

* ---- case xtcorr_re_effects
use "../data/sim_re_effects.dta", clear
xtset id t
generate double lx = L.x
scalar r_x_z_overall = .
scalar n_x_z_overall = .
preserve
quietly keep if !missing(x, z)
capture quietly correlate x z
if _rc == 0 scalar r_x_z_overall = r(rho)
if _rc == 0 scalar n_x_z_overall = r(N)
restore
scalar r_x_z_within = .
scalar n_x_z_within = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, fe clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_within = r(rho)
if _rc == 0 scalar n_x_z_within = r(N)
restore
scalar r_x_z_between = .
scalar n_x_z_between = .
preserve
quietly keep if !missing(x, z)
quietly xtdata x z, be clear
capture quietly correlate x z
if _rc == 0 scalar r_x_z_between = r(rho)
if _rc == 0 scalar n_x_z_between = r(N)
restore
scalar r_x_lx_overall = .
scalar n_x_lx_overall = .
preserve
quietly keep if !missing(x, lx)
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_overall = r(rho)
if _rc == 0 scalar n_x_lx_overall = r(N)
restore
scalar r_x_lx_within = .
scalar n_x_lx_within = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, fe clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_within = r(rho)
if _rc == 0 scalar n_x_lx_within = r(N)
restore
scalar r_x_lx_between = .
scalar n_x_lx_between = .
preserve
quietly keep if !missing(x, lx)
quietly xtdata x lx, be clear
capture quietly correlate x lx
if _rc == 0 scalar r_x_lx_between = r(rho)
if _rc == 0 scalar n_x_lx_between = r(N)
restore
scalar r_z_lx_overall = .
scalar n_z_lx_overall = .
preserve
quietly keep if !missing(z, lx)
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_overall = r(rho)
if _rc == 0 scalar n_z_lx_overall = r(N)
restore
scalar r_z_lx_within = .
scalar n_z_lx_within = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, fe clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_within = r(rho)
if _rc == 0 scalar n_z_lx_within = r(N)
restore
scalar r_z_lx_between = .
scalar n_z_lx_between = .
preserve
quietly keep if !missing(z, lx)
quietly xtdata z lx, be clear
capture quietly correlate z lx
if _rc == 0 scalar r_z_lx_between = r(rho)
if _rc == 0 scalar n_z_lx_between = r(N)
restore
mata: pyx_export("../out/xtcorr_re_effects__rho.csv", (st_numscalar("r_x_z_overall") \ st_numscalar("r_x_z_within") \ st_numscalar("r_x_z_between") \ st_numscalar("r_x_lx_overall") \ st_numscalar("r_x_lx_within") \ st_numscalar("r_x_lx_between") \ st_numscalar("r_z_lx_overall") \ st_numscalar("r_z_lx_within") \ st_numscalar("r_z_lx_between")))
mata: pyx_export("../out/xtcorr_re_effects__n.csv", (st_numscalar("n_x_z_overall") \ st_numscalar("n_x_z_within") \ st_numscalar("n_x_z_between") \ st_numscalar("n_x_lx_overall") \ st_numscalar("n_x_lx_within") \ st_numscalar("n_x_lx_between") \ st_numscalar("n_z_lx_overall") \ st_numscalar("n_z_lx_within") \ st_numscalar("n_z_lx_between")))

exit, clear
