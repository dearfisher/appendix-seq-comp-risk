# f and g functions for the alpha spending approach
# (See Jennison and Turnbull book Chapter 2013, p.9)

# rho-family spending functions (Kim-DeMets) for alpha and beta 
g <- function(I,rho,beta,Imax){
    beta*min(1,(I/Imax)^rho)
}
