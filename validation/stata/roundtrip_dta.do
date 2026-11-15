set more off
display "=== quarterly file written by pyxtpanel"
use "../out/rt_quarter.dta", clear
xtset
display "tdelta=" r(tdelta) " gaps=" r(gaps) " tmin=" r(tmins) " tmax=" r(tmaxs) " bal=" r(balanced)
display "=== delta(2) file"
use "../out/rt_delta.dta", clear
xtset
display "tdelta=" r(tdelta) " gaps=" r(gaps) " bal=" r(balanced)
list
