# This function finds the Cause Specific Hazard (CSH) ratio from :
# - the CIFs of main and competing events at time t among non treated
# - time t
# - the CIF of main at time t among non treated
# by assuming:
#  - constant CSH
#  - that the CSH of the competing event does not change with treatment.
#    as in  Pintilie Stat Med 2002 and Latouche 2007
# 
#

FindHRcsh2 <- function(CIF11,   # CIF main event, treated
                       CIF10,   # CIF main event, not treated
                       CIF20,   # CIF competing event, not treated
                       t,       # time at wich CIF are given
                       maxHR=100){ # maximum value for the HR, for the interval [0,maxHR] within wich we want to search the value of HR
    a10 <- - (CIF10*log(1 - CIF10 - CIF20)/(CIF10 + CIF20))/t
    a20 <- - (CIF20*log(1 - CIF10 - CIF20)/(CIF10 + CIF20))/t
    a21 <- a20
    a11 <- uniroot(function(x){CIF11 - (x/(x + a21))*(1 - exp(-(x + a21)*t))}, lower = 0, upper = maxHR*a10, tol = 1e-9)$root  
    CIF21 <- (a21/(a11 + a21))*(1 - exp(-(a11 + a21)*t))
    ## print(paste("To check",-CIF21*log(1 - CIF11 - CIF21)/(CIF11 + CIF21)/t,"and",a21))
    HR <- a11/a10
    res <- list(CSH10=a10,      # CSH main event control group
                CSH20=a20,      # CSH competing event control group
                CSH11=a11,      # CSH main event treatment group
                CSH21=a21,      # CSH competing event treatment group
                CIF10=CIF10,    # CIF main event at t for control group
                CIF20=CIF20,    # CIF competing event at t for control group
                CIF11=CIF11,    # CIF main event at t for treatment group
                CIF21=CIF21,    # CIF competing event at t for treatment group
                HR=HR           # Hazard ratio (for CSH) of main event i.e a11/a10
                )
    res
}
