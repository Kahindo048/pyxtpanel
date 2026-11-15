* L.3.1 / L.5.1 feasibility: which reference commands exist and what they return
foreach c in xtcd2 xtcsd xtserial xthst {
    capture which `c'
    display "`c' rc=" _rc
}
use "../data/sim_balanced.dta", clear
describe, short
xtset id t
xtdescribe, patterns(3)
quietly xtreg x z, fe
capture noisily xtcsd, pesaran
return list
capture noisily xtcd2 x z
return list
capture noisily xtcd2 x, noest
return list
capture noisily xtserial x z
return list
use "../data/sim_gaps.dta", clear
xtset id t
capture noisily xtserial x z
return list
capture noisily xtcd2 x, noest
return list
use "../data/sim_entry_exit.dta", clear
xtset id t
capture noisily xtserial x z
return list
capture noisily xtcd2 x, noest
return list
