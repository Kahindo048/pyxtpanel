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

* ---- case tsfill_gaps
use "../data/sim_gaps.dta", clear
xtset id t
tsfill
sort id t
mata: pyx_export("../out/tsfill_gaps__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_gaps__nrows.csv", st_nobs())

* ---- case tsfill_full_gaps
use "../data/sim_gaps.dta", clear
xtset id t
tsfill, full
sort id t
mata: pyx_export("../out/tsfill_full_gaps__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_full_gaps__nrows.csv", st_nobs())

* ---- case tsfill_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
tsfill
sort id t
mata: pyx_export("../out/tsfill_entry_exit__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_entry_exit__nrows.csv", st_nobs())

* ---- case tsfill_full_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
tsfill, full
sort id t
mata: pyx_export("../out/tsfill_full_entry_exit__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_full_entry_exit__nrows.csv", st_nobs())

* ---- case tsfill_weak
use "../data/sim_weak.dta", clear
xtset id t
tsfill
sort id t
mata: pyx_export("../out/tsfill_weak__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_weak__nrows.csv", st_nobs())

* ---- case tsfill_full_weak
use "../data/sim_weak.dta", clear
xtset id t
tsfill, full
sort id t
mata: pyx_export("../out/tsfill_full_weak__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_full_weak__nrows.csv", st_nobs())

* ---- case tsfill_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
tsfill
sort id t
mata: pyx_export("../out/tsfill_delta2__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_delta2__nrows.csv", st_nobs())

* ---- case tsfill_full_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
tsfill, full
sort id t
mata: pyx_export("../out/tsfill_full_delta2__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_full_delta2__nrows.csv", st_nobs())

* ---- case tsfill_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
tsfill
sort id t
mata: pyx_export("../out/tsfill_quarterly__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_quarterly__nrows.csv", st_nobs())

* ---- case tsfill_full_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
tsfill, full
sort id t
mata: pyx_export("../out/tsfill_full_quarterly__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_full_quarterly__nrows.csv", st_nobs())

* ---- case tsfill_monthly
use "../data/sim_monthly.dta", clear
xtset id t
tsfill
sort id t
mata: pyx_export("../out/tsfill_monthly__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_monthly__nrows.csv", st_nobs())

* ---- case tsfill_full_monthly
use "../data/sim_monthly.dta", clear
xtset id t
tsfill, full
sort id t
mata: pyx_export("../out/tsfill_full_monthly__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_full_monthly__nrows.csv", st_nobs())

* ---- case tsfill_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
tsfill
sort id t
mata: pyx_export("../out/tsfill_unsorted_ids__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_unsorted_ids__nrows.csv", st_nobs())

* ---- case tsfill_full_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
tsfill, full
sort id t
mata: pyx_export("../out/tsfill_full_unsorted_ids__data.csv", st_data(., ("id", "t", "x", "z")))
mata: pyx_export("../out/tsfill_full_unsorted_ids__nrows.csv", st_nobs())

exit, clear
