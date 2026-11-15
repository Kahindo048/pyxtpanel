* L.2.5 feasibility: xtunitroot fisher on unbalanced, short and gappy panels
foreach d in sim_entry_exit sim_gaps sim_sparse_t sim_one_obs_unit {
    use "../data/`d'.dta", clear
    quietly xtset id t
    display "=== `d'"
    quietly xtdescribe
    display "N=" r(N) " Tmin=" r(min) " Tmax=" r(max)
    capture noisily xtunitroot fisher x, dfuller lags(0)
    return list
}
use "../data/sim_ur_mixed.dta", clear
xtset id t
xtunitroot fisher y, dfuller lags(1) generate(pp)
return list
xtunitroot fisher y, dfuller lags(1) trend demean
xtunitroot fisher y, dfuller lags(2) drift
* per-unit dfuller for the first unit
dfuller y if id == 1, lags(1)
return list
dfuller y if id == 1, lags(1) trend
dfuller y if id == 1, lags(1) drift
