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

* ---- case xtsum_balanced_x
use "../data/sim_balanced.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_balanced_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_balanced_z
use "../data/sim_balanced.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_balanced_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_gaps_x
use "../data/sim_gaps.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_gaps_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_gaps_z
use "../data/sim_gaps.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_gaps_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_gaps_L_x
use "../data/sim_gaps.dta", clear
xtset id t
xtsum L.x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_gaps_L_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_gaps_D_z
use "../data/sim_gaps.dta", clear
xtset id t
xtsum D.z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_gaps_D_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_gaps_L2_z
use "../data/sim_gaps.dta", clear
xtset id t
xtsum L2.z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_gaps_L2_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_entry_exit_x
use "../data/sim_entry_exit.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_entry_exit_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_entry_exit_z
use "../data/sim_entry_exit.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_entry_exit_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_weak_x
use "../data/sim_weak.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_weak_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_weak_z
use "../data/sim_weak.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_weak_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_delta2_x
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_delta2_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_delta2_z
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_delta2_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_delta2_L_x
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
xtsum L.x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_delta2_L_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_quarterly_x
use "../data/sim_quarterly.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_quarterly_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_quarterly_z
use "../data/sim_quarterly.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_quarterly_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_quarterly_D_x
use "../data/sim_quarterly.dta", clear
xtset id t
xtsum D.x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_quarterly_D_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_unsorted_ids_x
use "../data/sim_unsorted_ids.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_unsorted_ids_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_unsorted_ids_z
use "../data/sim_unsorted_ids.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_unsorted_ids_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_one_obs_unit_x
use "../data/sim_one_obs_unit.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_one_obs_unit_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_one_obs_unit_z
use "../data/sim_one_obs_unit.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_one_obs_unit_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_single_unit_x
use "../data/sim_single_unit.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_single_unit_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_single_unit_z
use "../data/sim_single_unit.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_single_unit_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_single_period_x
use "../data/sim_single_period.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_single_period_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_single_period_z
use "../data/sim_single_period.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_single_period_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_constants_x
use "../data/sim_constants.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_constants_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_constants_cu
use "../data/sim_constants.dta", clear
xtset id t
xtsum cu
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_constants_cu__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_constants_ck
use "../data/sim_constants.dta", clear
xtset id t
xtsum ck
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_constants_ck__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_re_effects_x
use "../data/sim_re_effects.dta", clear
xtset id t
xtsum x
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_re_effects_x__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

* ---- case xtsum_re_effects_z
use "../data/sim_re_effects.dta", clear
xtset id t
xtsum z
scalar s_N = r(N)
scalar s_n = r(n)
scalar s_Tbar = r(Tbar)
scalar s_mean = r(mean)
scalar s_sd = r(sd)
scalar s_min = r(min)
scalar s_max = r(max)
scalar s_sd_b = r(sd_b)
scalar s_min_b = r(min_b)
scalar s_max_b = r(max_b)
scalar s_sd_w = r(sd_w)
scalar s_min_w = r(min_w)
scalar s_max_w = r(max_w)
mata: pyx_export("../out/xtsum_re_effects_z__r.csv", (st_numscalar("s_N") \ st_numscalar("s_n") \ st_numscalar("s_Tbar") \ st_numscalar("s_mean") \ st_numscalar("s_sd") \ st_numscalar("s_min") \ st_numscalar("s_max") \ st_numscalar("s_sd_b") \ st_numscalar("s_min_b") \ st_numscalar("s_max_b") \ st_numscalar("s_sd_w") \ st_numscalar("s_min_w") \ st_numscalar("s_max_w")))

exit, clear
