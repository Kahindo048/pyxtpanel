* L.6.2 feasibility: what Stata stores for Hansen / Sargan
clear all
set seed 12345
set obs 200
gen id = ceil(_n/10)
bysort id: gen t = _n
xtset id t
gen z1 = rnormal()
gen z2 = rnormal()
gen z3 = rnormal()
gen u = rnormal()
gen x = 0.6*z1 + 0.4*z2 + 0.2*z3 + 0.5*u + rnormal()
gen y = 1 + 0.8*x + u

display "=== ivregress 2sls + estat overid (Sargan)"
quietly ivregress 2sls y (x = z1 z2 z3)
estat overid
return list

display "=== ivregress gmm cluster + estat overid (Hansen J)"
quietly ivregress gmm y (x = z1 z2 z3), wmatrix(cluster id)
estat overid
return list

display "=== xtabond2 two-step robust"
gen yy = y
xtabond2 yy L.yy x, gmm(L.yy, lag(2 .)) iv(z1 z2) twostep robust nodiffsargan
ereturn list

display "=== is xtdpthresh installed?"
capture which xtdpthresh
display "rc = " _rc
