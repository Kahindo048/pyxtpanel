* L.6.3 feasibility: contents of e(diffsargan)
use "../data/sim_balanced.dta", clear
xtset id t
quietly xtabond2 x L.x z, gmm(x, lag(2 .)) iv(z) noleveleq noconstant twostep robust
matrix M = e(diffsargan)
matrix list M, format(%12.8f)
display "hansen full = " e(hansen) "  df " e(hansen_df)
use "../data/sim_gaps.dta", clear
xtset id t
quietly xtabond2 x L.x z, gmm(x, lag(2 .)) gmm(z, lag(1 2)) iv(z) noleveleq noconstant twostep robust
matrix M = e(diffsargan)
matrix list M, format(%12.8f)
display "groups: " e(diffgroup1) " | " e(diffgroup2)
