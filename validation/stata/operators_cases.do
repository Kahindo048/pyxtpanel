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

* ---- case operators_gaps
use "../data/sim_gaps.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_gaps__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_entry_exit__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_weak
use "../data/sim_weak.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_weak__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_delta2__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_quarterly__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_monthly
use "../data/sim_monthly.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_monthly__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_unsorted_ids__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_one_obs_unit__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_single_unit
use "../data/sim_single_unit.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_single_unit__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_single_period
use "../data/sim_single_period.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_single_period__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_all_missing_var
use "../data/sim_all_missing_var.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_all_missing_var__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

* ---- case operators_missing_keys
use "../data/sim_missing_keys.dta", clear
xtset id t
quietly generate double c1 = L.x
quietly generate double c2 = L2.x
quietly generate double c3 = L3.x
quietly generate double c4 = F.x
quietly generate double c5 = F2.x
quietly generate double c6 = D.x
quietly generate double c7 = D2.x
quietly generate double c8 = S2.x
quietly generate double c9 = S3.x
quietly generate double c10 = LD.x
quietly generate double c11 = L2D.x
quietly generate double c12 = FD.x
quietly generate double c13 = DS2.x
quietly generate double c14 = LS2.x
quietly generate double c15 = D2S2.x
quietly generate double c16 = L.z
quietly generate double c17 = F.z
quietly generate double c18 = D.z
quietly generate double c19 = L2.z
quietly generate double c20 = LD.z
quietly generate double c21 = S2.z
quietly generate double c22 = D2.z
drop if id >= . | t >= .
sort id t
mata: pyx_export("../out/operators_missing_keys__data.csv", st_data(., ("id", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")))

exit, clear
