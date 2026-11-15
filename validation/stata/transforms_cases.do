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

* ---- case transforms_fe_balanced
use "../data/sim_balanced.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_balanced__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_balanced
use "../data/sim_balanced.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_balanced__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_gaps
use "../data/sim_gaps.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_gaps__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_gaps
use "../data/sim_gaps.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_gaps__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_entry_exit__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_entry_exit__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_weak
use "../data/sim_weak.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_weak__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_weak
use "../data/sim_weak.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_weak__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_delta2__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_delta2__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_quarterly__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_quarterly__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_unsorted_ids__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_unsorted_ids__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_one_obs_unit__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_one_obs_unit__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_single_unit
use "../data/sim_single_unit.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_single_unit__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_single_unit
use "../data/sim_single_unit.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_single_unit__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_single_period
use "../data/sim_single_period.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_single_period__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_single_period
use "../data/sim_single_period.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_single_period__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_fe_re_effects
use "../data/sim_re_effects.dta", clear
xtset id t
xtdata x z, fe clear
sort id
mata: pyx_export("../out/transforms_fe_re_effects__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_be_re_effects
use "../data/sim_re_effects.dta", clear
xtset id t
xtdata x z, be clear
sort id
mata: pyx_export("../out/transforms_be_re_effects__data.csv", st_data(., ("id", "x", "z")))

* ---- case transforms_quasi
use "../data/sim_re_effects.dta", clear
xtset id t
quietly xtreg x z, re
scalar drift_u = abs(e(sigma_u) - 1.119347976677613)
scalar drift_e = abs(e(sigma_e) - 0.7748551012200685)
scalar su = 1.119347976677613
scalar se = 0.7748551012200685
bysort id: gen double Ti = _N
gen double theta = 1 - se / sqrt(se^2 + Ti * su^2)
quietly summarize theta
scalar drift_min = abs(r(min) - e(thta_min))
scalar drift_max = abs(r(max) - e(thta_max))
bysort id: egen double mx = mean(x)
bysort id: egen double mz = mean(z)
gen double xq = x - theta * mx
gen double zq = z - theta * mz
sort id t
mata: pyx_export("../out/transforms_quasi__data.csv", st_data(., ("id", "xq", "zq", "theta")))
mata: pyx_export("../out/transforms_quasi__sigma_drift.csv", (st_numscalar("drift_u") \ st_numscalar("drift_e") \ st_numscalar("drift_min") \ st_numscalar("drift_max")))

* ---- case transforms_twoway_balanced
use "../data/sim_balanced.dta", clear
xtset id t
quietly regress x i.id i.t
predict double ex, residuals
quietly regress z i.id i.t
predict double ez, residuals
sort id t
mata: pyx_export("../out/transforms_twoway_balanced__data.csv", st_data(., ("id", "ex", "ez")))

* ---- case transforms_twoway_gaps
use "../data/sim_gaps.dta", clear
xtset id t
quietly regress x i.id i.t
predict double ex, residuals
quietly regress z i.id i.t
predict double ez, residuals
sort id t
mata: pyx_export("../out/transforms_twoway_gaps__data.csv", st_data(., ("id", "ex", "ez")))

* ---- case transforms_twoway_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
quietly regress x i.id i.t
predict double ex, residuals
quietly regress z i.id i.t
predict double ez, residuals
sort id t
mata: pyx_export("../out/transforms_twoway_entry_exit__data.csv", st_data(., ("id", "ex", "ez")))

* ---- case transforms_twoway_weak
use "../data/sim_weak.dta", clear
xtset id t
quietly regress x i.id i.t
predict double ex, residuals
quietly regress z i.id i.t
predict double ez, residuals
sort id t
mata: pyx_export("../out/transforms_twoway_weak__data.csv", st_data(., ("id", "ex", "ez")))

* ---- case transforms_twoway_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
quietly regress x i.id i.t
predict double ex, residuals
quietly regress z i.id i.t
predict double ez, residuals
sort id t
mata: pyx_export("../out/transforms_twoway_delta2__data.csv", st_data(., ("id", "ex", "ez")))

* ---- case transforms_twoway_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
quietly regress x i.id i.t
predict double ex, residuals
quietly regress z i.id i.t
predict double ez, residuals
sort id t
mata: pyx_export("../out/transforms_twoway_quarterly__data.csv", st_data(., ("id", "ex", "ez")))

* ---- case transforms_twoway_re_effects
use "../data/sim_re_effects.dta", clear
xtset id t
quietly regress x i.id i.t
predict double ex, residuals
quietly regress z i.id i.t
predict double ez, residuals
sort id t
mata: pyx_export("../out/transforms_twoway_re_effects__data.csv", st_data(., ("id", "ex", "ez")))

* ---- case transforms_fod_balanced
use "../data/sim_balanced.dta", clear
xtset id t
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_balanced__data.csv", st_data(., ("id", "fod_x", "fod_z")))

* ---- case transforms_fod_gaps
use "../data/sim_gaps.dta", clear
xtset id t
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_gaps__data.csv", st_data(., ("id", "fod_x", "fod_z")))

* ---- case transforms_fod_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_entry_exit__data.csv", st_data(., ("id", "fod_x", "fod_z")))

* ---- case transforms_fod_weak
use "../data/sim_weak.dta", clear
xtset id t
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_weak__data.csv", st_data(., ("id", "fod_x", "fod_z")))

* ---- case transforms_fod_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_delta2__data.csv", st_data(., ("id", "fod_x", "fod_z")))

* ---- case transforms_fod_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_quarterly__data.csv", st_data(., ("id", "fod_x", "fod_z")))

* ---- case transforms_fod_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_unsorted_ids__data.csv", st_data(., ("id", "fod_x", "fod_z")))

* ---- case transforms_fod_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
gsort id -t
foreach v in x z {
    quietly generate byte ok_`v' = !missing(`v')
    quietly generate double v0_`v' = cond(ok_`v', `v', 0)
    quietly by id: generate double cs_`v' = sum(v0_`v')
    quietly by id: generate double cn_`v' = sum(ok_`v')
    quietly generate double ln_`v' = cn_`v' - ok_`v'
    quietly generate double fod_`v' = sqrt(ln_`v' / (ln_`v' + 1)) * (`v' - (cs_`v' - v0_`v') / ln_`v') if ok_`v' & ln_`v' > 0
}
sort id t
mata: pyx_export("../out/transforms_fod_one_obs_unit__data.csv", st_data(., ("id", "fod_x", "fod_z")))

exit, clear
