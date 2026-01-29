
use "D:\Temp\MCP\DID\did-toolbox\Data\quota_stata.dta", clear
cd D:\Temp\MCP\DID\did-toolbox\examples\
cap log close 
log using quota_logStata, replace

encode country, generate(country_code)
tsset country_code year
xtreg womparl quota, fe vce(cluster country_code)


* no covariates & dependant: quota dummy
sdid womparl country year quota, vce(bootstrap) seed(1213)

drop if lngdp==.
sdid womparl country year quota, vce(bootstrap) seed(1213) covariates(lngdp, projected)


* no covariates & dependant: lnmmrt
tsset country_code year
xtreg lnmmrt quota if lnmmrt~=., fe vce(cluster country_code)
sdid lnmmrt country year quota if lnmmrt~=., vce(bootstrap) seed(1213)


log  close  
print quota_logStata.smcl