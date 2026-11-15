use "../data/invest.dta", clear
xtset firm year
xtdescribe
display "balanced=" r(balanced) " N=" r(N) " T=" r(T_max)
