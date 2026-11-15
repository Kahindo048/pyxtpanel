* L.3.2-3.4 feasibility
foreach c in xttest2 xtcsd xtcd2 {
    capture which `c'
    display "`c' rc=" _rc
}
use "../data/sim_csd.dta", clear
xtset id t
quietly xtreg y x, fe
display "g_min=" e(g_min) " N_g=" e(N_g) " T=" 12
xtcsd, friedman
return list
xtcsd, frees
return list
xtcsd, pesaran abs show
return list
capture noisily xttest2
return list
use "../data/sim_csd_unbalanced.dta", clear
xtset id t
quietly xtreg y x, fe
display "g_min=" e(g_min)
capture noisily xtcsd, friedman
return list
capture noisily xtcsd, frees
return list
