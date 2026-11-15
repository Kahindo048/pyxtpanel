set more off
* ---- panel de test : desequilibre, trou interne (unite 2 saute t=3), unite 3 courte
clear
input id t x z
1 1 1 10
1 2 2 11
1 3 3 12
1 4 4 13
2 1 5 .
2 2 5 15
2 4 7 17
3 2 2 20
3 3 2 21
end
display "=== A. xtset + return list"
xtset id t
return list
display "=== B. xtdescribe + return list"
xtdescribe
return list
display "=== C. xtsum x z + return list (x)"
xtsum x
return list
display "T-bar etc"
quietly xtsum z
display "z sd_w=" %18.15g r(sd_w) " N=" r(N) " n=" r(n)
* sd within a la main : sqrt( sum (x - xbar_i)^2 / (N-1) )
preserve
keep if z<.
bysort id: egen m = mean(z)
gen d2 = (z-m)^2
quietly summ d2
display "manual sdw (N-1)=" %18.15g sqrt(r(sum)/(r(N)-1)) "  (N-n)=" %18.15g sqrt(r(sum)/(r(N)-3))
restore
display "=== D. xttab x + r(results)"
xttab x
matrix list r(results), format(%12.6f)
return list
display "=== E. xttrans x (gap du temps pontee ?)"
xttrans x, freq
return list
display "=== F. tsfill"
preserve
tsfill
list, sepby(id)
restore
preserve
tsfill, full
list, sepby(id)
restore
display "=== G. xtdata fe / be"
preserve
xtdata x z, fe clear
list
restore
preserve
xtdata x z, be clear
list
restore
