# This function finds the CSH for the main and competing event from :
#  - CIF of main event at t for control group
#  - CIF competing event at t for control group
#  - time t
#  - Hazard ratio (for CSH) of main event 
# by assuming:
#  - constant CSH
#  - that the CSH of the competing event does not change with treatment.
#
# The function also compute CIF at t for the treatment group
#
# This follows more or less the approach explained in
# Schulgen et al Contemporary Clinical Trials 2005
  
CIFsandHRtoAllCSH <- function(CIF10,  # CIF main event at t for control group
                              CIF20,   # CIF competing event at t for control group
                              t,     # time t
                              HR     # Hazard ratio (for CSH) of main event 
                              ){
    # compute everything
    a10 <- - (CIF10*log(1 - CIF10 - CIF20)/(CIF10 + CIF20))/t
    a20 <- - (CIF20*log(1 - CIF10 - CIF20)/(CIF10 + CIF20))/t
    a11 <- a10*HR
    a21 <- a20
    CIF11 <- (a11/(a11 + a21))*(1 - exp(-(a11 + a21)*t))
    CIF21 <- (a21/(a11 + a21))*(1 - exp(-(a11 + a21)*t))
    # save results
    res <- list(a10=a10,       # CSH main event control group
                a20=a20,       # CSH competing event control group
                a11=a11,       # CSH main event treatment group
                a21=a21,       # CSH competing event treatment group
                CIF10=CIF10,   # CIF main event at t for control group
                CIF20=CIF20,   # CIF competing event at t for control group
                CIF11=CIF11,   # CIF main event at t for treatment group
                CIF21=CIF21,   # CIF competing event at t for treatment group
                HR=HR          # Hazard ratio (for CSH) of main event i.e a11/a10
                )
    res
}
