/*analysisDenuncias.do v1.00     damiancclarke             yyyy-mm-dd:2020-08-25
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This code generates summary statistics, two-way FE models, and panel event study
analyses of the impact of quarantine imposition on "Calls to number 149 of 
Carabineros (FonoFamilia) classified as VIF, with communal and monthly breakdown 
from 2019 to date". It requires the following ados installed from the SSC.
 
 > eventdd
 > did_multiplegt
 > estout

All globals are set in section 0 to import data and export results and log files.
Replicating results should simply require altering these globals to addresses on 
your own computer, and installing the above mentioned ados.
*/

vers 15
clear all
set more off
cap log close

set matsize 1000
*-------------------------------------------------------------------------------
*--- (0) Globals, set-up
*-------------------------------------------------------------------------------
global ROOT "C:/Users/danie/Desktop/Proyectos/Asistente/replication"
global ROOT "/home/damian/investigacion/2020/IPV_COVID19/replication"

global DAT "$ROOT/data/VIF/dta"
global OUT "$ROOT/results/vif"
global LOG "$ROOT/log"


cap mkdir "$OUT"
cap mkdir "$LOG"
log using "$LOG/analysisVIF.txt", text replace
foreach folder in eventdd areg did_multiplegt descriptives {
    cap mkdir "$OUT/`folder'"
}
foreach ado in eventdd did_multiplegt estout {
    cap which `ado'
    if _rc!=0 ssc install `ado'
}

*-------------------------------------------------------------------------------
*--- (1) Open data and creating new variables
*-------------------------------------------------------------------------------
use "$DAT/VIF.dta", clear

local outcomes VIF VIF1 VIF2 VIF3

local outcomespc
foreach var of varlist `outcomes' {
    gen `var'pc = `var'/population*100000
    local outcomespc `outcomespc' `var'pc
}

gen postCOVID = year == 2020 & month>=3
areg VIFpc postCOVID, abs(comuna)
areg VIFpc postCOVID quarantine, abs(comuna)


bys comuna (t): gen n = _n
bys comuna (t): egen minn = min(n) if quarantine==1
bys comuna: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ

levelsof comuna, local(comunas)
foreach tp in externo interno {
    foreach c of local comunas {
        if "`c'"!="11302"&"`c'"!="12202" {
            qui sum mobility_`tp' if t==721&comuna==`c'
            replace mobility_`tp'=`r(mean)' if t<721&comuna==`c'
        }
    }
}
lab var quarantine "Quarantine / Calls (arbitrary scale)"
lab var VIFpc "Quarantine / Calls (arbitrary scale)"
lab var mobility_interno "Mobility"
lab var mobility_externo "Mobility"
replace VIFpc = VIFpc/40

#delimit ;
twoway  line quarantine t if nom_com=="Ñuñoa", lwidth(thick) yaxis(2)
|| connected mobility_interno t if  nom_com=="Ñuñoa", ylabel(4(1)8) ms(Dh)
|| connected mobility_externo t if  nom_com=="Ñuñoa", scheme(plottig) ms(Sh)
||      line VIFpc t if nom_com=="Ñuñoa", yaxis(2) lpattern(dash) lwidth(medthick)
legend(order(1 "Quarantine (left axis)" 2 "Calls to Helpline (left axis)" 3 "Internal Mobility (right axis)"
             4 "External Mobility (right axis)") position(6) row(2))
xtitle("Month") xline(722.5, lcolor(red));
graph export "$OUT/descriptives/generalSetup_wCalls.eps", replace;

twoway  line quarantine t if nom_com=="Ñuñoa", lwidth(thick) yaxis(2)
|| connected mobility_interno t if  nom_com=="Ñuñoa", ylabel(4(1)8) ms(Dh)
|| connected mobility_externo t if  nom_com=="Ñuñoa", scheme(plottig) ms(Sh) 
legend(order(1 "Quarantine (left axis)" 2 "Internal Mobility (right axis)"
             3 "External Mobility (right axis)") position(6) row(1))
xtitle("Month");
#delimit cr
graph export "$OUT/descriptives/generalSetup.eps", replace


*-------------------------------------------------------------------------------
*--- (1b) Descriptives
*-------------------------------------------------------------------------------
gen qmonth = month if quarantine == 1
bys comuna: egen minQ = min(qmonth)
bys comuna: egen maxQ = max(qmonth)

gen qend = maxQ+12 if maxQ!=9
gen timeToExit = n - qend

set scheme plottig
#delimit ;
eventdd VIFpc  i.t i.comuna if qstart!=., timevar(timeToExit) ci(rcap) lags(18) 
leads(6) baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black)) 
graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6" -3
                "-3" 0 "0" 2 "2" 4 "4" 6 "{&ge} 6")
         scheme(plottig) xtitle("Months Relative to Quarantine Exit")
         ytitle("Calls to #149 per 100,000 people"))
endpoints_op(ms(Dh) mc(midblue));
graph export "$OUT/eventdd/quarantineExit.eps", replace;
#delimit cr

foreach num of numlist 2(1)17 {
    gen lagQ`num' = timeToQ==-`num'
    gen lagE`num' = timeToExit==-`num'
}
gen lagQ18 = timeToQ<=-18
gen lagE18 = timeToExit<=-18
foreach num of numlist 0(1)5 {
    gen leadQ`num' = timeToQ==`num'
    gen leadE`num' = timeToExit==`num'
}
gen leadQ6 = timeToQ>=6 & timeToQ!=.
gen leadE6 = timeToExit>=6 & timeToExit!=.

areg VIFpc i.t lagQ* leadQ* lagE* leadE*, cluster(comuna) abs(comuna)


preserve



gen postQuarantine = month>maxQ&year==2020
reghdfe VIFpc quarantine postQuarantine mobility_*, absorb(t comuna) cluster(comuna)
reghdfe VIFpc quarantine , absorb(t comuna) cluster(comuna)
reghdfe VIFpc mobility_* , absorb(t comuna) cluster(comuna)
reghdfe VIFpc mobility_* quarantine, absorb(t comuna) cluster(comuna)




collapse (sum) VIF population, by(month year maxQ)
gen VIFpc = VIF/population*100000
bys maxQ (year month): gen time = _n
set scheme plotplainblind
#delimit ;
twoway connected VIFpc time if maxQ==4
||     connected VIFpc time if maxQ==5
||     connected VIFpc time if maxQ==7
||     connected VIFpc time if maxQ==8
||     connected VIFpc time if maxQ==.,
legend(order(1 "Out Apr" 2 "Out May" 3 "Out Jul"
             4 "Never Out" 5 "Never In"))
xlabel(1 "Jan 2019" 4 "Apr 2019" 7 "Jul 2019" 10 "Oct 2019" 13 "Jan 2020"
       16 "Apr 2020" 18 "Jun 2020" 20 "Aug 2020", angle(45)) xtitle("")
xline(15, lcolor(red)) ytitle("Calls to #149 for DV per 100,000 Inhabitants");
graph export "$OUT/descriptives/in_out.eps", replace;
#delimit cr
restore


preserve
collapse (sum) VIF1 VIF2 VIF3, by(t)
gen time = _n
gen total = VIF1+VIF2+VIF3
tab time
count
drop if time>21
rename VIF1 economic
rename VIF2 physical
rename VIF3 psychological
set scheme plotplainblind
#delimit ;
twoway connected total time, lwidth(thick)
||  connected economic time, lwidth(medthick)
||  connected physical time, lwidth(medthick)
||  connected psychological time, lwidth(medthick)
xlabel(1 "Jan 2019" 4 "Apr 2019" 7 "Jul 2019" 10 "Oct 2019" 13 "Jan 2020"
       16 "Apr 2020" 19 "Jul 2020" 21 "Sep 2020", angle(45)) xtitle("")
ytitle("Calls to #149 for DV")
legend(order(1 "Total" 2 "Economic" 3 "Physical" 4 "Psychological"))
xline(15, lcolor(red));
graph export "$OUT/descriptives/calls149.eps", replace;
#delimit cr
restore

preserve
lab var VIFpc "All Calls to \#149 per 100,000 Inhabitants"
lab var VIF1pc "Calls to \#149 for Economic Violence per 100,000"
lab var VIF2pc "Calls to \#149 for Physical Violence per 100,000"
lab var VIF3pc "Calls to \#149 for Psychological Violence per 100,000"

#delimit ;
estpost sum VIFpc VIF1pc VIF2pc VIF3pc;
estout using "$OUT/descriptives/Summary149.tex", replace label style(tex)
cells("count mean(fmt(2)) sd(fmt(2)) min(fmt(1)) max(fmt(1))")
collabels(, none) mlabels(, none);
#delimit cr
restore

*-------------------------------------------------------------------------------
*--- (2) Event Study
*-------------------------------------------------------------------------------
local fes i.t i.comuna 
generate upper=15

foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt
        local gn 
    }
    if "`wt'"=="yes" {
        local opt [aw=population]
        local gn _Wt 
    }
    local j = 1
    foreach var of varlist `outcomespc' {
        if `j'==1 local x = -4
        if `j'==2 local x = -0.2
        if `j'==3 local x = -2
        if `j'==4 local x = -1
        
        #delimit ;
        eventdd `var' `fes' `opt' if timeToExit<0|timeToExit==., timevar(timeToQ) ci(rcap) lags(18) 
        leads(6) baseline(-7) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6" -3
                        "-3" 0 "0" 2 "2" 4 "4" 6 "{&ge} 6")
                 scheme(plottig) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people")
                 xline(-7) text(`x' -4 "Partially Post-COVID", size(small))
                 text(`x' 2 "Post-Quarantine", size(small)))
        endpoints_op(ms(Dh) mc(midblue));
        graph export "$OUT/eventdd/event1`gn'_`var'.eps", replace;
	
        eventdd `var' `fes' population `opt', timevar(timeToQ) ci(rcap)
        lags(18) leads(6) baseline(-7) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6" -3
                        "-3" 0 "0" 1 "+1" 2 "+2" 3 "3" 4 "4" 5 "5" 6 "6+")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people")
                 xline(-7) text(`x' -4 "Partially Post-COVID", size(small))
                 text(`x' 2 "Post-Quarantine", size(small)));
        graph export "$OUT/eventdd/event2`gn'_`var'.eps", replace;
	
        eventdd `var' `fes' mobility_externo mobility_interno `opt',
        timevar(timeToQ) ci(rcap) lags(18) leads(6)
        baseline(-7) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6" -3
                        "-3" 0 "0" 1 "+1" 2 "+2" 3 "3" 4 "4" 5 "5" 6 "6+")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people")
                 xline(-7) text(`x' -4 "Partially Post-COVID", size(small))
                 text(`x' 2 "Post-Quarantine", size(small)));
        graph export "$OUT/eventdd/event3`gn'_`var'.eps", replace;
	
        eventdd `var' `fes' population mobility_externo mobility_interno `opt'
        if timeToExit<0|timeToExit==.,
        timevar(timeToQ) ci(rcap) lags(18) leads(6)
        baseline(-7) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6" -3
                        "-3" 0 "0" 1 "+1" 2 "+2" 3 "3" 4 "4" 5 "5" 6 "6+")
                 scheme(plottig) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people")
                 xline(-7) text(`x' -4 "Partially Post-COVID", size(small))
                 text(`x' 2 "Post-Quarantine", size(small)));
        graph export "$OUT/eventdd/event4`gn'_`var'.eps", replace;
        #delimit cr
        local ++j
    }
    graph drop _all
}

#delimit ;
eventdd VIFpc `fes' if qstart==15|qstart==16|qstart==., timevar(timeToQ) ci(rcap)
lags(15) leads(5) baseline(-2) coef_op(ms(Dh)) accum ci_op(lcolor(black))
graph_op(xlabel(-15 "{&le} -15" -12 "-12" -9 "-9" -6 "-6" -3
                "-3" 0 "0" 1 "+1" 2 "+2" 3 "3" 4 "4" 5 "5" 6 "6+")
         scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
         ytitle("Calls to #149 per 100,000 people"));
graph export "$OUT/eventdd/event_VIF_earlyQs.eps", replace;

gen schoolClosed=t>722;
eventdd VIFpc `fes' mobility_*, timevar(timeToQ) ci(rcap)
lags(15) leads(6) baseline(-7) coef_op(ms(Dh)) accum ci_op(lcolor(black))
graph_op(xlabel(-15 "{&le} -15" -12 "-12" -9 "-9" -6 "-6" -3
                "-3" 0 "0" 1 "+1" 2 "+2" 3 "3" 4 "4" 5 "5" 6 "6+")
         scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
         ytitle("Calls to #149 per 100,000 people"));
graph export "$OUT/eventdd/event_VIF_schoolclosed.eps", replace;
#delimit cr
exit


*-------------------------------------------------------------------------------
*--- (3) Two way FEs
*-------------------------------------------------------------------------------
xtset comuna t
local fes i.t i.comuna 
local se abs(comuna) cluster(comuna) 

lab var quarantine "Quarantine Imposed"
lab var mobility_interno "Internal Movement in Municipality"
lab var mobility_externo "External Movement in Municipality"
foreach var of varlist `outcomespc' {
    *Sin Pesos de Poblacion
    eststo: areg `var' `fes' quarantine, `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' population quarantine, `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_externo mobility_interno quarantine, `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_externo mobility_interno population quarantine, `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    *Con Pesos de Poblacion
    eststo: areg `var' `fes' quarantine [aw=population], `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' population quarantine [aw=population], `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_externo mobility_interno quarantine [aw=population], `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_externo mobility_interno population quarantine [aw=population], `se'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    local ests est1 est2 est3 est4 est5 est6 est7 est8
    #delimit ;
    esttab `ests' using "$OUT/areg/DD_`var'.tex",
    b(%-9.3f) se(%-9.3f) noobs nonotes nogaps
    keep(mobility_externo mobility_interno population quarantine)
    mlabels(, none) nonumbers style(tex) fragment replace noline label
    stats(N mean, fmt(%9.0gc %5.3f)
          label("\\ Observations" "Mean of Dependent Variable"))
    starlevel ("*" 0.10 "**" 0.05 "***" 0.01);
    #delimit cr
    estimates clear
}

*-------------------------------------------------------------------------------
*--- (4) Sharp Difference-in-Difference
*-------------------------------------------------------------------------------
local fes comuna t
local eopts placebo(5) dynamic(3) breps(50) cluster(comuna)
local eopts3 placebo(3) dynamic(3) breps(50) cluster(comuna)

foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt
        local gn 
    }
    if "`wt'"=="yes" {
        local opt weight(population)
        local gn _Wt 
    }
    foreach var of varlist `outcomespc' {
        did_multiplegt `var' `fes' quarantine, `eopts' `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did1_quar`gn'_`var'.eps", replace
        
        did_multiplegt `var' `fes' quarantine, `eopts' controls(population) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did2_quar`gn'_`var'.eps", replace
        
        did_multiplegt `var' `fes' quarantine, `eopts3' controls(mobility_externo mobility_interno) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did3_quar`gn'_`var'.eps", replace
        
        did_multiplegt `var' `fes' quarantine, `eopts3' controls(mobility_externo mobility_interno population) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did4_quar`gn'_`var'.eps", replace
    }
    graph drop _all
}


log close
