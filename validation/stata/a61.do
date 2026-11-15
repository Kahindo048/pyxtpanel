* L.6.1 feasibility: AR(m) in xtabond2 and estat abond
use "../data/sim_balanced.dta", clear
xtset id t
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan twostep robust
display "xtabond2 ar1=" e(ar1) " p=" e(ar1p) "  ar2=" e(ar2) " p=" e(ar2p)
quietly xtabond x, lags(1) noconstant twostep vce(robust)
estat abond, artests(2)
return list
quietly xtabond2 x L.x, gmm(x, lag(2 .)) noleveleq noconstant nodiffsargan robust
display "xtabond2 onestep ar1=" e(ar1) " ar2=" e(ar2)
quietly xtabond x, lags(1) noconstant vce(robust)
estat abond, artests(2)
return list
