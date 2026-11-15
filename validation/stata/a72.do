use "../data/sim_one_obs_unit.dta", clear
xtset id t
egen g = group(id)
gen byte okeq = !missing(D.x)
bysort id: egen neq = total(okeq)
tab neq
levelsof id if g == 29, local(i)
display "unit with g=29: id `i'"
list id t x if g == 29
capture noisily dfuller x if id == `i', lags(0) regress
display "rc = " _rc
