# f and g functions for the alpha spending approach
# (See Jennison and Turnbull book Chapter 2013, p.9)

# rho-family spending functions (Kim-DeMets) for alpha and beta 
f <- function(I,rho,alpha,Imax){
  alpha*min(1,(I/Imax)^rho)
}
