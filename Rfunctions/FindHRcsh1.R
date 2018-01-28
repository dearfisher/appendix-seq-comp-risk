# This function finds the Cause Specific Hazard (CSH) ratio from :
# - the cumulative incidence functions of main
#   and competing at time t among treated and non treated
# - the time t
# This follows an approach explained in Pintilie Stat Med 2002,
# which essentially assumes
# - constant baseline
# - proportional hazards
# for all CSH.
  
FindHRcsh1 <- function(CIF11,   # CIF main event, treated
                       CIF21,   # CIF competing event, treated
                       CIF10,   # CIF main event, not treated
                       CIF20,   # CIF competing event, not treated
                       t){      # time at wich CIF are given
    a11 <- - (log(1 - CIF11 - CIF21)/(1 + CIF21/CIF11))/t
    ## a21 <- a11*CIF21/CIF11
    a10 <- - (log(1 - CIF10 - CIF20)/(1 + CIF20/CIF10))/t
    ## a20 <- a11*CIF20/CIF10
    HR <- a11/a10
    HR
}
