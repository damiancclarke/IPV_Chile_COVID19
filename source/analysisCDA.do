/*analysisDenuncias.do v1.00     damiancclarke             yyyy-mm-dd:2020-08-31
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This code generates summary statistics, two-way FE models, and panel event study
analyses of the impact of quarantine imposition on "Casa de Acogida". It req-
uires the following ados installed from the SSC.

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
set matsize 3000

*-------------------------------------------------------------------------------
*--- (0) Globals, set-up
*-------------------------------------------------------------------------------
global ROOT "C:/Users/danie/Desktop/Proyectos/Asistente/replication"
global ROOT "/home/damian/investigacion/2020/IPV_COVID19/replication"

global DAT "$ROOT/data/CDA/dta"
global OUT "$ROOT/results/cda"
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
use "$DAT/CDA.dta", clear

drop if t==21996 //21mar
drop if t==21995 //22mar

gen ingresspc = (ingreso_m+ingreso_n)/population*100000
gen occupancypc  = (ingreso_m/(ingreso_m+cupos_m))/population*100000

local outcomespc ingresspc occupancypc

bys Region (t): gen n = _n
bys Region (t): egen minn = min(n) if quarantine==1
bys Region: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ

levelsof Region, local(region)
foreach tp in externo interno {
    foreach r of local region {
        qui sum mobility_`tp' if t>=21971&t<=21989&Region==`r'
        replace mobility_`tp'=`r(mean)' if t<21971&Region==`r'
    }
}

*-------------------------------------------------------------------------------
*--- (2) Event Study
*-------------------------------------------------------------------------------
*By Day
local fes i.t i.Region
              
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
	if "`var'"=="ingresspc"   local et="Ingress of Residents per 100,000 people"
	if "`var'"=="occupancypc" local et="Occupancy of Residents per 100,000 people"
        
        #delimit ;
        eventdd `var' `fes' `opt', timevar(timeToQ) ci(rcap) 
		lags(50) leads(50) wboot baseline(-1) coef_op(ms(Dh)) 
		wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-50 "{&le} -50" -25 "-25" 0 "0" 25 "25" 50 "{&ge} 50")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event1`gn'_`var'.eps", replace;
    	
        eventdd `var' `fes' population `opt', timevar(timeToQ) ci(rcap) 
		lags(50) leads(50) wboot baseline(-1) coef_op(ms(Dh)) 
		wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-50 "{&le} -50" -25 "-25" 0 "0" 25 "25" 50 "{&ge} 50")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event2`gn'_`var'.eps", replace;
    
    	eventdd `var' `fes' mobility_ext mobility_int `opt',
        timevar(timeToQ) ci(rcap) lags(50) leads(50) wboot baseline(-1) 
        coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-50 "{&le} -50" -25 "-25" 0 "0" 25 "25" 50 "{&ge} 50")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event3`gn'_`var'.eps", replace;
    	
    	eventdd `var' `fes' population mobility_ext mobility_int `opt', 
        timevar(timeToQ) ci(rcap) lags(50) leads(50) wboot baseline(-1) 
        coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-50 "{&le} -50" -25 "-25" 0 "0" 25 "25" 50 "{&ge} 50")
                 scheme(s1mono) xtitle("Days Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event4`gn'_`var'.eps", replace;
    	#delimit cr
    }
    graph drop _all
}

*By Week
preserve
drop timeToQ
sort Region t
by Region: gen time = _n
gen week = ceil(time/7)
collapse population mobility_ext mobility_int quarantine (sum) `outcomespc', by(Region week)
gen minn = week if quarantine!=0
bys Region: egen qstart = min(minn)
gen timeToQ = week-qstart
local fes i.week i.Region

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
	if "`var'"=="ingresspc"   local et="Ingress of Residents per 100,000 people"
	if "`var'"=="occupancypc" local et="Occupancy of Residents per 100,000 people"
        
        #delimit ;
        eventdd `var' `fes' `opt', timevar(timeToQ) ci(rcap) 
		lags(10) leads(9) wboot baseline(-1) coef_op(ms(Dh)) 
		wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-10 "{&le} -10" -5 "-5" 0 "0" 5 "5" 9 "{&ge} 9")
                 scheme(s1mono) xtitle("Weeks Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event1`gn'_`var'_wk.eps", replace;
    	
        eventdd `var' `fes' population `opt', timevar(timeToQ) ci(rcap) 
		lags(10) leads(9) wboot baseline(-1) coef_op(ms(Dh)) 
		wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-10 "{&le} -10" -5 "-5" 0 "0" 5 "5" 9 "{&ge} 9")
                 scheme(s1mono) xtitle("Weeks Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event2`gn'_`var'_wk.eps", replace;
    
    	eventdd `var' `fes' mobility_ext mobility_int `opt',
        timevar(timeToQ) ci(rcap) lags(10) leads(9) wboot baseline(-1) 
        coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-10 "{&le} -10" -5 "-5" 0 "0" 5 "5" 9 "{&ge} 9")
                 scheme(s1mono) xtitle("Weeks Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event3`gn'_`var'_wk.eps", replace;
    	
    	eventdd `var' `fes' population mobility_ext mobility_int `opt', 
        timevar(timeToQ) ci(rcap) lags(10) leads(9) wboot baseline(-1) 
        coef_op(ms(Dh)) wboot_op(seed(1213)) accum ci_op(lcolor(black))
        graph_op(xlabel(-10 "{&le} -10" -5 "-5" 0 "0" 5 "5" 9 "{&ge} 9")
                 scheme(s1mono) xtitle("Weeks Relative to Quarantine Imposition")
                 ytitle("`et'"));
        graph export "$OUT/eventdd/event4`gn'_`var'_wk.eps", replace;
    	#delimit cr
    }
    graph drop _all
}
restore

*-------------------------------------------------------------------------------
*--- (3) Two way FEs
*-------------------------------------------------------------------------------
**DANIEL, ¿Aquí hay alguna manera de automatizar la inclusión de los intervalos
** de confianza de boottest en la tabla de regresión?

xtset Region t
local fes i.t i.Region
local wt [aw=population]
local se abs(Region) cluster(Region)

foreach q of varlist quarantine PropPopQuar {
    foreach var of varlist `outcomespc' {
        *Sin Pesos de Poblacion
        eststo: areg `var' `fes' `q', `se' 
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
        eststo: areg `var' `fes' population `q', `se'
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
        
        eststo: areg `var' `fes' mobility_ext mobility_int `q', `se'
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
        eststo: areg `var' `fes' mobility_ext mobility_int population `q', `se'
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
        *Con Pesos de Poblacion
        eststo: areg `var' `fes' `q' `wt', `se'
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
        eststo: areg `var' `fes' population `q' `wt', `se'
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
        eststo: areg `var' `fes' mobility_ext mobility_int `q' `wt', `se'
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
	
        eststo: areg `var' `fes' mobility_ext mobility_int population `q' `wt', `se'
        sum `var' if e(sample)==1
        estadd scalar mean=r(mean)
        boottest `q'
    
        local ests est1 est2 est3 est4 est5 est6 est7 est8
        #delimit ;
        esttab `ests' using "$OUT/areg/DD_`q'_`var'.tex",
        b(%-9.3f) se(%-9.3f) noobs nonotes nogaps
        keep(mobility_externo mobility_interno population `q') 
        mlabels(, none) nonumbers style(tex) fragment replace noline label
        starlevel ("*" 0.10 "**" 0.05 "***" 0.01);
        #delimit cr
        estimates clear
    }
    graph drop _all
}

*-------------------------------------------------------------------------------
*--- (4) Sharp Difference-in-Difference
*-------------------------------------------------------------------------------
sort Region t
by Region: gen time = _n
gen week = ceil(time/7)
collapse population mobility_ext mobility_int quarantine PropPopQuar (sum) ingresspc occupancypc, by(Region week)
replace quarantine=1 if quarantine!=0
local fes Region week

foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt placebo(5) dynamic(3) breps(50) cluster(Region) 
        local gn 
    }
    if "`wt'"=="yes" {
        local opt placebo(5) dynamic(3) breps(50) cluster(Region) weight(population)
        local gn _Wt 
    }
    foreach tv in quarantine PropPopQuar {
        if "`tv'"=="quarantine"  local v _quar
        if "`tv'"=="PropPopQuar" local v _pop
        foreach var of varlist `outcomespc' {
            did_multiplegt `var' `fes' `tv', `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did1`v'`gn'_`var'.eps", replace
            
            did_multiplegt `var' `fes' `tv', controls(population) `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did2`v'`gn'_`var'.eps", replace
            
            did_multiplegt `var' `fes' `tv', controls(mobility_ext mobility_int) `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did3`v'`gn'_`var'.eps", replace
            
            did_multiplegt `var' `fes' `tv', controls(mobility_ext mobility_int population) `opt'
            ereturn list
            graph export "$OUT/did_multiplegt/did4`v'`gn'_`var'.eps", replace
        }
    }
    graph drop _all
}

exit


*--------------
use "$DAT/CDA.dta", clear

drop if t==21996 //21mar
drop if t==21995 //22mar

gen ingresspc = (ingreso_m+ingreso_n)/population*100000
gen occupancypc  = (ingreso_m/(ingreso_m+cupos_m))/population*100000

local outcomespc ingresspc occupancypc

bys Region (t): gen n = _n
bys Region (t): egen minn = min(n) if quarantine==1
bys Region: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ

levelsof Region, local(region)
foreach tp in externo interno {
    foreach r of local region {
        qui sum mobility_`tp' if t>=21971&t<=21989&Region==`r'
        replace mobility_`tp'=`r(mean)' if t<21971&Region==`r'
    }
}

areg ingresspc i.t i.Region quarantine, abs(Region) cluster(Region)
boottest quarantine, nogr
*matrix def a=r(CI)
 
eststo est1, addscalars(ci1 r(CI)): areg ingresspc i.t i.Region quarantine, abs(Region) cluster(Region) 

#delimit ;
esttab est1, b(%-9.3f) se(%-9.3f) noobs nonotes nogaps
        keep(quarantine) 
        mlabels(, none) nonumbers style(tex) fragment replace noline label
        starlevel ("*" 0.10 "**" 0.05 "***" 0.01) scalars(ci1);
#delimit cr



