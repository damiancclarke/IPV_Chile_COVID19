/*analysisDenuncias.do v1.00     damiancclarke             yyyy-mm-dd:2020-09-03
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This code generates summary statistics, two-way FE models, and panel event study
analyses of the impact of quarantine imposition on "Complaints made to Carabineros 
from 2018, with communal and daily breakdown". It requires the following ados 
installed from the SSC.
 
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
*set matsize ....

*-------------------------------------------------------------------------------
*--- (0) Globals, set-up
*-------------------------------------------------------------------------------
global ROOT "C:/Users/danie/Desktop/Proyectos/Asistente/replication/data"
*global ROOT "/home/damian/investigacion/2020/IPV_COVID19/replication"

global DAT "$ROOT/VIF_Delitos/dta"
global OUT "$ROOT/Results/delitos"
global LOG "$ROOT/Results/log"

log using "$LOG/analysisDelitos.txt", text replace
cap mkdir $OUT

*-------------------------------------------------------------------------------
*--- (1) Open data and creating new variables
*-------------------------------------------------------------------------------
use "$DAT/Denuncias.dta", clear

local outcomes vmamujer vmahombre

local outcomespc
foreach var of varlist `outcomes' {
    gen `var'pc = `var'/population*100000
	local outcomespc `outcomespc' `var'pc
}

bys comuna (t): gen n = _n
bys comuna (t): egen minn = min(n) if quarantine==1
bys comuna: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ

drop if comuna==.|t==. //eliminar missing

levelsof comuna, local(comunas)
foreach tp in externo interno{
    foreach c of local comunas{
	    if "`c'"!="11302"&"`c'"!="12202"{
	        sum mobility_`tp' if t>=721&t<=725&comuna==`c'
	        replace mobility_`tp'=`r(mean)' if t<721&comuna==`c'
	    }
    }
}

*-------------------------------------------------------------------------------
*--- (2) Event Study 
*-------------------------------------------------------------------------------
*Quarantine is the event
local fes i.t i.comuna 

foreach wt in no yes {
    if "`wt'"=="no"  {
    local opt
    local gn 
    }
	if "`wt'"=="yes" {
    local opt [aw=population]
    local gn _Wt 
    }
foreach var of varlist `outcomespc' {
    #delimit ;
    eventdd `var' `fes' `opt', 
	timevar(timeToQ) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event1`gn'_`var'.eps", replace;
	
    eventdd `var' `fes' population `opt', 
	timevar(timeToQ) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event2`gn'_`var'.eps", replace;

	eventdd `var' `fes' mobility_ext mobility_int `opt', 
	timevar(timeToQ) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event3`gn'_`var'.eps", replace;
	
	eventdd `var' `fes' population mobility_ext mobility_int `opt', 
	timevar(timeToQ) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2" 3 "+3")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event4`gn'_`var'.eps", replace;
	#delimit cr
}
    graph drop _all
}

*Exit quarantine is the event
bys comuna: egen maxQ = max(n) if quarantine==1
bys comuna: egen qend = min(maxQ)
gen timeToExit = n-qend
tab timeToExit

local fes i.t i.comuna 

foreach wt in no yes {
    if "`wt'"=="no"  {
    local opt
    local gn 
    }
	if "`wt'"=="yes" {
    local opt [aw=population]
    local gn _Wt 
    }
foreach var of varlist `outcomespc' {
    #delimit ;
    eventdd `var' `fes' `opt', 
	timevar(timeToExit) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event1`gn'_`var'_exit.eps", replace;
	
    eventdd `var' `fes' population `opt', 
	timevar(timeToExit) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event2`gn'_`var'_exit.eps", replace;

	eventdd `var' `fes' mobility_ext mobility_int `opt', 
	timevar(timeToExit) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event3`gn'_`var'_exit.eps", replace;
	
	eventdd `var' `fes' population mobility_ext mobility_int `opt', 
	timevar(timeToExit) ci(rcap) lags(170) leads(3)
    baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
    graph_op(xlabel(-170 "{&le} -170" -130 "-130" -90 "-90" -40 "-40" -10 "-10" 0 "0" 1 "+1" 2 "+2")
             scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
             ytitle("Crime Rate (per capita)"));
    graph export "$OUT/eventdd/event4`gn'_`var'_exit.eps", replace;
	#delimit cr
}
    graph drop _all
}

*-------------------------------------------------------------------------------
*--- (3) Two way FEs
*-------------------------------------------------------------------------------
xtset comuna t
local fes i.t i.comuna 
local wt [aw=population]
local ops abs(comuna) cluster(comuna) 


foreach var of varlist `outcomespc' {
    dis "here"
    *Sin Pesos de Poblacion
    eststo: areg `var' `fes' quarantine, `ops' 
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
	
	eststo: areg `var' `fes' population quarantine, `ops'
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    boottest `q'
	
	eststo: areg `var' `fes' mobility_ext mobility_int quarantine, `ops'
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
	
	eststo: areg `var' `fes' mobility_ext mobility_int population quarantine, `ops'
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
	
    *Con Pesos de Poblacion
	eststo: areg `var' `fes' quarantine `wt', `ops'
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
	
	eststo: areg `var' `fes' population quarantine `wt', `ops'
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
	
	eststo: areg `var' `fes' mobility_ext mobility_int quarantine `wt', `ops'
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
	
	eststo: areg `var' `fes' mobility_ext mobility_int population quarantine `wt', `ops'
	sum `var' if e(sample)==1
    estadd scalar mean=r(mean)

    local ests est1 est2 est3 est4 est5 est6 est7 est8
	#delimit ;
	esttab `ests' using "$OUT/areg/DD_`var'.tex",
	b(%-9.3f) se(%-9.3f) noobs keep(mobility_ext mobility_int population quarantine) 
	nonotes nogaps mlabels(, none) nonumbers style(tex) fragment replace noline label
	starlevel ("*" 0.10 "**" 0.05 "***" 0.01);
	#delimit cr
	estimates clear
}

*-------------------------------------------------------------------------------
*--- (4) Sharp Difference-in-Difference
*-------------------------------------------------------------------------------
local fes comuna t
local eopts1 placebo(5) dynamic(3) placebo(5) cluster(comuna)

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
        did_multiplegt `var' `fes' quarantine, `eopts1' `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did1_quar`gn'_`var'.eps", replace

        did_multiplegt `var' `fes' quarantine, `eopts1' `opt' controls(population)
        ereturn list
        graph export "$OUT/did_multiplegt/did2_quar`gn'_`var'.eps", replace
        
        did_multiplegt `var' `fes' quarantine, `eopts1' `opt' controls(mobility_ext mobility_int) 
        ereturn list
        graph export "$OUT/did_multiplegt/did3_quar`gn'_`var'.eps", replace
        
        did_multiplegt `var' `fes' quarantine, `eopts1' `opt' controls(mobility_ext mobility_int population) 
        ereturn list
        graph export "$OUT/did_multiplegt/did4_quar`gn'_`var'.eps", replace
    }
    graph drop _all
}
exit 
