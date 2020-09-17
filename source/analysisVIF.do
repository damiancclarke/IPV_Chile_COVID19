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

*-------------------------------------------------------------------------------
*--- (0) Globals, set-up
*-------------------------------------------------------------------------------
global ROOT "C:/Users/danie/Desktop/Proyectos/Asistente/replication"

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

*-------------------------------------------------------------------------------
*--- (2) Event Study
*-------------------------------------------------------------------------------
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
        eventdd `var' `fes' `opt', timevar(timeToQ) ci(rcap) lags(18) 
		leads(3) baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6" -3
                        "-3" 0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people"));
        graph export "$OUT/eventdd/event1`gn'_`var'.eps", replace;
	
	eventdd `var' `fes' population `opt', timevar(timeToQ) ci(rcap)
        lags(18) leads(3) baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6"
                        -3 "-3" 0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people"));
        graph export "$OUT/eventdd/event2`gn'_`var'.eps", replace;
	
	eventdd `var' `fes' mobility_externo mobility_interno `opt',
        timevar(timeToQ) ci(rcap) lags(18) leads(3)
        baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6"
                        -3 "-3" 0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people"));
        graph export "$OUT/eventdd/event3`gn'_`var'.eps", replace;
	
	eventdd `var' `fes' population mobility_externo mobility_interno `opt',
        timevar(timeToQ) ci(rcap) lags(18) leads(3)
        baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-18 "{&le} -18" -15 "-15" -12 "-12" -9 "-9" -6 "-6" -3
                        "-3" 0 "0" 1 "+1" 2 "+2" 3 "+3")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Calls to #149 per 100,000 people"));
        graph export "$OUT/eventdd/event4`gn'_`var'.eps", replace;
        #delimit cr
    }
    graph drop _all
}

*-------------------------------------------------------------------------------
*--- (3) Two way FEs
*-------------------------------------------------------------------------------
xtset comuna t
local fes i.t i.comuna 
local se abs(comuna) cluster(comuna) 

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
    b(%-9.3f) se(%-9.3f) noobs keep(mobility_externo mobility_interno population quarantine) nonotes nogaps
    mlabels(, none) nonumbers style(tex) fragment replace noline label
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
