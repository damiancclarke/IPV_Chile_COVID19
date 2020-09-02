/*analysisDenuncias.do v1.00     damiancclarke             yyyy-mm-dd:2020-08-26
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
set matsize 2000


**TODO: CHANGE DESCRIPTIVE GRAPH LABELS TO ENGLISH
*-------------------------------------------------------------------------------
*--- (0) Globals, set-up
*-------------------------------------------------------------------------------
global ROOT "/home/damian/investigacion/2020/IPV_COVID19/replication/"

global DAT "$ROOT/data/SPD/dta"
global OUT "$ROOT/results/spd"
global LOG "$ROOT/log"



cap mkdir "$OUT"
cap mkdir "$LOG"
log using "$LOG/analysisSPD.txt", text replace
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
use "$DAT/SPD.dta", clear

local outcomes denuncias victimas

local outcomespc
foreach var of varlist `outcomes' {
    gen `var'pc = `var'/poblacion*100000
	local outcomespc `outcomespc' `var'pc
}

bys comuna (t): gen n = _n
bys comuna (t): egen minn = min(n) if quarantine==1
bys comuna: egen qstart = min(minn)
gen timeToQ = n-qstart
tab timeToQ

***DANIEL: mejor re-emplazamos con promedios en una linea base.  Por ejemplo,
***re-emplazar cada comuna con el valor promedio durante el periodo pre-covid
***(especificamente me refiero al periodo pre-covid donde tenemos datos de movilidad...)
replace mobility_externo=0 if mobility_externo==. //REVISAR
replace mobility_interno=0 if mobility_interno==. //REVISAR
drop if comuna==.|t==. //eliminar datos que no entran en el analisis

*-------------------------------------------------------------------------------
*--- (2) Descriptives
*-------------------------------------------------------------------------------
preserve
collapse (sum) denuncias, by(month year)
rename month orden
reshape wide denuncias, i(orden) j(year)

graph bar (sum) denuncias2019 denuncias2020, ///
blabel(total, position(outside) format(%5.0f)) ///
over(orden, relabel(1 "Enero" 2 "Febrero" 3 "Marzo" 4 "Abril" 5 "Mayo")) ///
bargap(5) bar(1, color(ebblue*1))  bar(2, colo(red*1.5)) ///
legend(label(1 "2019") label(2 "2020") ///
       pos(12) col(2)  size(vsmall) region(lcolor(white))) ///
ylabel(0(1000)10000, ///
       labsize(vsmall) angle(horizontal) glcolor(black*0.4) glpattern(solid) glwidth (vvthin) valuelabel format(%5.0f)) ///
ytitle("Cantidad de Denunicas (Violencia Intrafamiliar)") ///
graphregion(fcolor(white) color(black))
graph export "$OUT/descriptives/Denuncias_month_20192020.eps", replace
restore


preserve
collapse (sum) denuncias victimas_f victimas_m poblacion, by(month year t)
replace denuncias = denuncias/poblacion*10000
replace victimas_f = victimas_f/poblacion*10000
replace victimas_m = victimas_m/poblacion*10000
tsset t
replace t = t-365 if year>2018
replace t = t-365 if year>2019

#delimit ;
tsline denuncias if year==2018, lcolor(gs13) lpattern(solid) ||
tsline denuncias if year==2019, lcolor(gs13) lpattern(solid) ||
tsline denuncias if year==2020, lcolor(red) lpattern(solid) lwidth(medthick)
legend(order(1 "Año 2018" 2 "Año 2019" 3 "Año 2020"))
ytitle("Denuncias por 10,000 Habitantes") xtitle("")
xlabel(21185 "1 Enero" 21216 "1 Febrero" 21244 "1 Marzo"
       21275 "1 Abril" 21305 "1 Mayo" 21336 "1 Junio", angle(45));
graph export "$OUT/descriptives/denunciasNacional.eps", replace;

drop if t==.;
tsline victimas_f if year==2018, lcolor(gs13) lpattern(solid) ||
tsline victimas_f if year==2019, lcolor(gs13) lpattern(solid) ||
tsline victimas_f if year==2020, lcolor(gs3) lpattern(solid)
legend(order(1 "Año 2018" 2 "Año 2019" 3 "Año 2020"))
ytitle("Victimas Femininas por 10,000 Habitantes") xtitle("")
xlabel(21185 "1 Enero" 21216 "1 Febrero" 21244 "1 Marzo"
       21275 "1 Abril" 21305 "1 Mayo" 21336 "1 Junio", angle(45));
graph export "$OUT/descriptives/victimasFNacional.eps", replace;

tsline victimas_m if year==2018, lcolor(gs13) lpattern(solid) ||
tsline victimas_m if year==2019, lcolor(gs13) lpattern(solid) ||
tsline victimas_m if year==2020, lcolor(gs3) lpattern(solid)
legend(order(1 "Año 2018" 2 "Año 2019" 3 "Año 2020"))
ytitle("Victimas Masculinas por 10,000 Habitantes") xtitle("")
xlabel(21185 "1 Enero" 21216 "1 Febrero" 21244 "1 Marzo"
       21275 "1 Abril" 21305 "1 Mayo" 21336 "1 Junio", angle(45));
graph export "$OUT/descriptives/victimasMNacional.eps", replace;
#delimit cr
restore

preserve
collapse (sum) denuncias victimas_f victimas_m poblacion, by(month year t region)
replace denuncias = denuncias/poblacion*10000
replace victimas_f = victimas_f/poblacion*10000
replace victimas_m = victimas_m/poblacion*10000


#delimit ;
label define R 1 "Región I" 2 "Región II" 3 "Región III" 4 "Región IV"
5 "Región V" 6 "Región VI" 7 "Región VII" 8 "Región VIII"
9 "Región IX" 10 "Región X" 11 "Región XI" 12 "Región XII"
13 "Región XIII" 14 "Región XIV" 15 "Región XV" 16 "Región XVI";
#delimit cr
label val region R

drop if denuncias>0.7
#delimit ;
twoway connected denuncias t if year==2020, by(region, note("")) c(l)
m(none) xtitle("") xlabel(,angle(45)) xline(21990) lwidth(medthick)
ytitle("Denuncias por 10,000 Habitantes");
graph export "$OUT/descriptives/denunciasRegion149.eps", replace;

drop if t==.;
twoway connected victimas_f t if year==2020, by(region, note("")) c(l)
m(none) xtitle("") xlabel(,angle(45)) xline(21990) lwidth(medthick)
ytitle("Denuncias por 10,000 Habitantes");
graph export "$OUT/descriptives/victimasRegion149.eps", replace;

twoway connected victimas_m t if year==2020, by(region, note("")) c(l)
m(none) xtitle("") xlabel(,angle(45)) xline(21990) lwidth(medthick)
ytitle("Denuncias por 10,000 Habitantes");
graph export "$OUT/descriptives/victimasMRegion149.eps", replace;
#delimit cr
restore

***DANIEL: Si no es muy complicado, podemos agregar los graficos comunales
***que compara tasas de 2019 y 2020 y que tiene la linea de 45 grados que
***tenemos en el archivo SPDdesc.do.  Avísame si no tienes el archivo y te
***lo envio!  


*-------------------------------------------------------------------------------
*--- (3) Event Study
*-------------------------------------------------------------------------------
***DANIEL: En este caso como el event study es por día, mejor hacemos un análisis
***con más lags y leads.  Va a continuación un codigo que tengo con un análisis
***diario y semanal.  ¿Podemos implementar ambos aquí.  Puede ser en dos bloques
***de code distinto si es más simple.
/*
drop if comuna==""
keep if year==2020
gen rate = denuncias/poblacionINE*10000
sort comuna month day
by comuna: gen time = _n
gen date = time if quarantine==1
bys comuna: egen adopDate = min(date)
gen timeToAdop = time-adopDate
egen cid = group(comuna)

graph set eps fontface "Times New Roman"
eventdd rate i.cid i.time [aw=poblacionINE], timevar(timeToAdop) ci(rcap) cluster(cid) lags(50) leads(50) graph_op(legend(off) scheme(s1mono) xtitle("Días a Cuarentena") ytitle("Denuncias por 10,000 Habitantes") xlabel(-50 "{&le} -50" -25 "-25" 0 "0" 25 "25" 50 "{&ge} 50")) coef_op(ms(Oh) mc(black)) ci_op(color(gs12)) endpoints_op(ms(Oh) mc(black)) accum
graph export eventDenuncias_baseline.eps, replace


drop time

sort comuna month day
by comuna: gen time = _n
gen week = ceil(time/7)
collapse poblacionINE quarantine (sum) denuncias, by(comuna week)
gen date = week if quarantine!=0
bys comuna: egen adopDate = min(date)
gen timeToAdop = week-adopDate
egen cid = group(comuna)
gen rate = denuncias/poblacionINE*10000

eventdd rate i.cid i.week [aw=poblacionINE], timevar(timeToAdop) ci(rcap) cluster(cid) lags(10) leads(7) graph_op(legend(off) scheme(s1mono) xtitle("Semanas a Cuarentena") ytitle("Denuncias por 10,000 Habitantes") xlabel(-10 "{&le} -10" -5 "-5" 0 "0" 5 "5" 7 "{&ge} 7")) coef_op(ms(Oh) mc(black)) ci_op(color(gs12)) endpoints_op(ms(Oh) mc(black)) accum
graph export eventDenuncias_baseline_weeks.eps, replace

*/

foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt
        local gn 
    }
    if "`wt'"=="yes" {
        local opt [aw=poblacion]
        local gn _Wt 
    }
    foreach var of varlist `outcomespc' {
        #delimit ;
        eventdd `var' i.t i.comuna `opt', timevar(timeToQ) ci(rcap) lags(20)
        leads(3) baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1"
                        0 "0" 1 "+1" 2 "+2" 3 "+3" 5 "+5")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Crime Rate (per capita)"));
        graph export "$OUT/eventdd/event1`gn'_`var'.eps", replace;
	
        eventdd `var' i.t i.comuna poblacion `opt', timevar(timeToQ) ci(rcap)
        lags(20) leads(3) baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1"
                        0 "0" 1 "+1" 2 "+2" 3 "+3" 5 "+5")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Crime Rate (per capita)"));
        graph export "$OUT/eventdd/event2`gn'_`var'.eps", replace;
        
	eventdd `var' i.t i.comuna mobility_externo mobility_interno `opt',
        timevar(timeToQ) ci(rcap) lags(20) leads(3)
        baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1"
                        0 "0" 1 "+1" 2 "+2" 3 "+3" 5 "+5")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Crime Rate (per capita)"));
        graph export "$OUT/eventdd/event3`gn'_`var'.eps", replace;
	
	eventdd `var' i.t i.comuna poblacion mobility_externo mobility_interno `opt',
        timevar(timeToQ) ci(rcap) lags(20) leads(3)
        baseline(-1) coef_op(ms(Dh)) accum ci_op(lcolor(black))
        graph_op(xlabel(-20 "{&le} -20" -15 "-15" -10 "-10" -5 "-5" -1 "-1" 0
                        "0" 1 "+1" 2 "+2" 3 "+3" 5 "+5")
                 scheme(s1mono) xtitle("Months Relative to Quarantine Imposition")
                 ytitle("Crime Rate (per capita)"));
        graph export "$OUT/eventdd/event4`gn'_`var'.eps", replace;
	#delimit cr
    }
    graph drop _all
}

*-------------------------------------------------------------------------------
*--- (4) Two way FEs
*-------------------------------------------------------------------------------
xtset comuna t
local ops abs(comuna) cluster(comuna) 
local fes i.t i.comuna
local wt [aw=poblacion]

foreach var of varlist `outcomespc' {
    dis "here"
    *Sin Pesos de Poblacion
    eststo: areg `var' `fes' quarantine, `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' poblacion quarantine, `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_ext mobility_int quarantine, `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_ext mobility_int poblacion quarantine, `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    *Con Pesos de Poblacion
    eststo: areg `var' `fes' quarantine `wt', `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' poblacion quarantine `wt', `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_ext mobility_int quarantine `wt', `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)
    
    eststo: areg `var' `fes' mobility_ext mobility_int poblacion quarantine `wt', `ops'
    sum `var' if e(sample)==1
    estadd scalar mean=r(mean)

    local ests est1 est2 est3 est4 est5 est6 est7 est8
    #delimit ;
    esttab `ests' using "$OUT/areg/DD_`q'_`var'.tex",
    b(%-9.3f) se(%-9.3f) noobs keep(`q') nonotes nogaps mlabels(, none) 
    nonumbers style(tex) fragment replace noline label
    starlevel ("*" 0.10 "**" 0.05 "***" 0.01);
    #delimit cr
    estimates clear
}

*-------------------------------------------------------------------------------
*--- (5) Sharp Difference-in-Difference
*-------------------------------------------------------------------------------
***DANIEL: Aquí propongo que trabajemos a nivel de semana...  Si realizamos algo
*** como el proceso de hacer un estudio de ventos semanal, nos va a dar mucho
*** mas estabilidad en los 'lags' y 'leads'.  Si no es claro, me avises! GRACIAS!

local fes comuna t
foreach wt in no yes {
    if "`wt'"=="no"  {
        local opt placebo(5) dynamic(3) breps(50) cluster(comuna)
        local gn 
    }
    if "`wt'"=="yes" {
        local opt placebo(5) dynamic(3) breps(50) cluster(comuna) weight(poblacion)
        local gn _Wt 
    }
    foreach var of varlist `outcomespc' {
        did_multiplegt `var' `fes' quarantine, `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did1_quar`gn'_`var'.eps", replace

        did_multiplegt `var' `fes' quarantine, controls(poblacion) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did2_quar`gn'_`var'.eps", replace
        
        did_multiplegt `var' `fes' quarantine, controls(mobility_ext mobility_int) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did3_quar`gn'_`var'.eps", replace
        
        did_multiplegt `var' `fes' quarantine, controls(mobility_ext mobility_int poblacion) `opt'
        ereturn list
        graph export "$OUT/did_multiplegt/did4_quar`gn'_`var'.eps", replace
    }
    graph drop _all
}




 
 
 
 
