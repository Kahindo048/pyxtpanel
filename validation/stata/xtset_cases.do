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

* ---- case xtset_balanced
use "../data/sim_balanced.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_balanced__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_gaps
use "../data/sim_gaps.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_gaps__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_entry_exit
use "../data/sim_entry_exit.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_entry_exit__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_weak
use "../data/sim_weak.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_weak__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_delta2
use "../data/sim_delta2.dta", clear
xtset id t, delta(2)
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_delta2__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_quarterly
use "../data/sim_quarterly.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_quarterly__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_monthly
use "../data/sim_monthly.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_monthly__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_unsorted_ids
use "../data/sim_unsorted_ids.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_unsorted_ids__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_constants
use "../data/sim_constants.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_constants__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_one_obs_unit
use "../data/sim_one_obs_unit.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_one_obs_unit__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_single_unit
use "../data/sim_single_unit.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_single_unit__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_single_period
use "../data/sim_single_period.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_single_period__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_all_missing_var
use "../data/sim_all_missing_var.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_all_missing_var__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_missing_keys
use "../data/sim_missing_keys.dta", clear
xtset id t
scalar tmin = r(tmin)
scalar tmax = r(tmax)
scalar tdelta = r(tdelta)
scalar gaps = r(gaps)
scalar imin = r(imin)
scalar imax = r(imax)
scalar bal = cond("`r(balanced)'" == "strongly balanced", 1, cond("`r(balanced)'" == "weakly balanced", 2, 3))
scalar nobs = _N
preserve
drop if id >= . | t >= .
scalar nobs = _N
bysort id: gen byte _first = (_n == 1)
quietly summarize _first
scalar nunits = r(sum)
bysort t: gen byte _firstt = (_n == 1)
quietly summarize _firstt
scalar nperiods = r(sum)
restore
mata: pyx_export("../out/xtset_missing_keys__xtset.csv", (st_numscalar("tmin") \ st_numscalar("tmax") \ st_numscalar("tdelta") \ st_numscalar("gaps") \ st_numscalar("imin") \ st_numscalar("imax") \ st_numscalar("bal") \ st_numscalar("nobs") \ st_numscalar("nunits") \ st_numscalar("nperiods")))

* ---- case xtset_refused_duplicates
use "../data/sim_duplicates.dta", clear
capture xtset id t
scalar rc = _rc
mata: pyx_export("../out/xtset_refused_duplicates__rc.csv", (st_numscalar("rc")))

* ---- case xtset_refused_non_integer_time
use "../data/sim_non_integer_time.dta", clear
capture xtset id t
scalar rc = _rc
mata: pyx_export("../out/xtset_refused_non_integer_time__rc.csv", (st_numscalar("rc")))

exit, clear
