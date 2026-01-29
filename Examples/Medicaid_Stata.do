clear
cd "d:\Temp\MCP\DID\did-toolbox\Examples\"
cap log close
log using medicaid, replace
* just 2013 vs 2014
import excel "D:\Projekte\basic_econometrics\Diff_in_Diff\DiDToolbox\Medicaid_own.xlsx", sheet("Sheet1") firstrow

tsset county_code year_code
keep if yaca==.|(yaca==2014|yaca>2019)

xtreg crude_rate_20_64 treatPost  i.year if year_code>=2013 & year_code<=2014 , fe vce(cluster county_code)

* all years
import excel "D:\Projekte\basic_econometrics\Diff_in_Diff\DiDToolbox\Medicaid_own.xlsx", sheet("Sheet1") firstrow clear
tsset county_code year_code


xtdidregress (crude_rate_20_64) (treatPost), group(county_code) time(year_code)
estat bdecomp

** Wooldridge
xthdidregress twfe (crude_rate_20_64) (treatPost), group(county_code) 
estat aggregation, cohort 
estat aggregation, time


** Callaway/Sant'Anna DR
xthdidregress twfe (crude_rate) (treatPost), group(county_code) 
estat aggregation, cohort
estat aggregation, time

xthdidregress aipw (crude_rate) (treatPost perc_female perc_white perc_hispanic unemp_rate_pc poverty_rate median_income_k), group(county_code) 
estat aggregation, cohort
estat aggregation, time


** CH?
gen Dhat = treatPost
replace Dhat=. if Dhat==0
*did_imputation crude_rate_20_64 county_code year_code Dhat, controls(perc_female perc_white perc_hispanic unemp_rate_pc poverty_rate median_income_k) *autosample

log close

print medicaid.smcl