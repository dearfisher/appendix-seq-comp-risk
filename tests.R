rm(list=ls())

## {{{ Technicalities
## {{{ load packages
library(mvtnorm)
library(survival)
library(cmprsk)
library(Publish)
library(gsDesign)
library(testthat)
## }}}

## {{{ source functions
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}

# path to where Rfunctions have been copied
sourceDir("path/to/Rfunctions/")
## }}}

## {{{ Create .Random.seed
runif(1) # .Random.seed is created when you first call a random number generator. 
## }}}


## }}}


##----------


## {{{ Test values of Table 7.6 (p. 164) Jennison and Turnbull's book 
test_that("Table 7.6 p164 : rho=2, K=4, beta=0.2, alpha=0.05",{
    MyResults <- planboundaries(rho=2,  
                                  alpha=0.05,
                                  beta=0.2,
                                  K=4,         
                                  theta=log(2),
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=FALSE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )
    expect_equal(MyResults$coef,1.087,tolerance=5e-4)
})

test_that("Table 7.6 p. 164 : rho=2, K=5, beta=0.2, alpha=0.05",{
    MyResults <- planboundaries(rho=2,  
                                  alpha=0.05,
                                  beta=0.2,
                                  K=5,         
                                  theta=log(2),
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=FALSE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )
    expect_equal(MyResults$coef,1.098,tolerance=5e-4)
})


test_that("Table 7.6 p. 164 : rho=3, K=5, beta=0.2, alpha=0.05",{
    MyResults <- planboundaries(rho=3,  
                              alpha=0.05,
                              beta=0.2,
                              K=5,         
                              theta=log(2),
                              abseps = 1e-08,
                              binding=TRUE,
                              direction="smaller",
                              Trace=FALSE,
                              nWhileMax=30,
                              toldiff=1e-06,
                              mycoefMax= 1.5,       
                              mycoefL=1             
                              )
    expect_equal(MyResults$coef,1.045,tolerance=5e-4)
})

test_that("Table 7.6 p. 164 : rho=2, K=3, beta=0.2, alpha=0.05",{
MyResults <- planboundaries(rho=2,  
                              alpha=0.05,
                              beta=0.2,
                              K=3,         
                              theta=log(2),
                              abseps = 1e-08,
                              binding=TRUE,
                              direction="smaller",
                              Trace=FALSE,
                              nWhileMax=30,
                              toldiff=1e-06,
                              mycoefMax= 1.2,       
                              mycoefL=1             
                              )
    expect_equal(MyResults$coef,1.070,tolerance=5e-4)
})


test_that("Table 7.6 p. 164 : rho=2, K=2, beta=0.2, alpha=0.05",{

    MyResults <- planboundaries(rho=2,  
                                  alpha=0.05,
                                  beta=0.20,
                                  K=2,         
                                  theta=log(2),
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=FALSE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )
    expect_equal(MyResults$coef,1.043,tolerance=5e-4)
})


test_that("Table 7.6 p. 164 : rho=2, K=5, beta=0.1, alpha=0.10",{
    MyResults <- planboundaries(rho=2,  
                                  alpha=0.05,
                                  beta=0.10,
                                  K=5,         
                                  theta=1,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=TRUE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )    
    expect_equal(MyResults$coef,1.100,tolerance=5e-4)
})

test_that("Table 7.6 p. 164 : rho=2, K=5, beta=0.1, alpha=0.05",{
    MyResults <- planboundaries(rho=2,  
                                  alpha=0.05,
                                  beta=0.10,
                                  K=5,         
                                  theta=1,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=TRUE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )    
    expect_equal(MyResults$coef,1.100,tolerance=5e-4)
})


test_that("Table 7.6 p. 164 : rho=3, K=5, beta=0.05, alpha=0.05",{
    MyResults <- planboundaries(rho=3,  
                                  alpha=0.05,
                                  beta=0.05,
                                  K=5,         
                                  theta=3,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=TRUE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )    
    expect_equal(MyResults$coef,1.050,tolerance=5e-4)
})


test_that("Table 7.6 p. 164 : rho=3, K=6, beta=0.05, alpha=0.05",{
    MyResults <- planboundaries(rho=3,  
                                  alpha=0.05,
                                  beta=0.05,
                                  K=6,         
                                  theta=1,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=TRUE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )    
    expect_equal(MyResults$coef,1.055,tolerance=5e-4)
})


test_that("Table 7.6 p. 164 : rho=3, K=10, beta=0.05, alpha=0.05",{
    MyResults <- planboundaries(rho=3,  
                                  alpha=0.05,
                                  beta=0.05,
                                  K=10,         
                                  theta=1,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=TRUE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.5,       
                                  mycoefL=1             
                                  )    
    expect_equal(MyResults$coef,1.069,tolerance=5e-4)
})

test_that("Table 7.6 p. 164 : rho=3, K=12, beta=0.10, alpha=0.05",{
    MyResults <- planboundaries(rho=3,  
                                  alpha=0.05,
                                  beta=0.10,
                                  K=12,         
                                  theta=0.5,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=TRUE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.2,       
                                  mycoefL=1             
                                  )    
    expect_equal(MyResults$coef,1.070,tolerance=5e-4)
})

test_that("Table 7.6 p. 164 : rho=2, K=9, beta=0.20, alpha=0.05",{
    MyResults <- planboundaries(rho=2,  
                                alpha=0.05,
                                beta=0.20,
                                K=9,         
                                theta=0.25,
                                abseps = 1e-08,
                                binding=TRUE,
                                direction="smaller",
                                Trace=TRUE,
                                nWhileMax=30,
                                toldiff=1e-06,
                                mycoefMax= 1.2,       
                                mycoefL=1             
                                )    
    expect_equal(MyResults$coef,1.120,tolerance=5e-4)
})
## }}}


## {{{ Test reproduce Table 25.2 of Chapter Jennison and Turnbull 2013


test_that("Reproduce Table 25.2 of Chapter Jennison and Turnbull 2013 (1)",{
    MyResults <- planboundaries(rho=2,  
                                  alpha=0.025,
                                  beta=0.2,
                                  K=5,         
                                  theta=0.5,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=FALSE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.2,       
                                  mycoefL=1             
                                  )
    ## print(MyResults)
    expect_equal(MyResults$Imax,34.48,tolerance=5e-3)
})

test_that("Reproduce Table 25.2 of Chapter Jennison and Turnbull 2013 (2)",{
    MyResults <- planboundaries(rho=2,  
                                  alpha=0.025,
                                  beta=0.2,
                                  K=5,         
                                  theta=0.5,
                                  abseps = 1e-08,
                                  binding=FALSE,
                                  direction="smaller",
                                  Trace=FALSE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.2,       
                                  mycoefL=1             
                                  )
    expect_equal(MyResults$Imax,35.58,tolerance=5e-3)
})


## }}}


## {{{ Test reproduce Table 25.3 of Chapter Jennison and Turnbull 2013

test_that("Test reproduce Table 25.3 of Chapter Jennison and Turnbull 2013",{
    Planned <- planboundaries(rho=2,  
                              alpha=0.025,
                              beta=0.2,
                              K=5,         
                              theta=0.5,
                              abseps = 1e-08,
                              binding=TRUE,
                              direction="smaller",
                              Trace=FALSE,
                              nWhileMax=30,
                              toldiff=1e-06,
                              mycoefMax= 1.2,       
                              mycoefL=1             
                              )

    boundaries1 <- InterimAnalysis(SeqCR.object=Planned,
                                   hatbeta.k=-1.04*sqrt(5.43),
                                   hatse.k=1/sqrt(5.43))

    boundaries2 <- InterimAnalysis(SeqCR.object=boundaries1,
                                   hatbeta.k=-1.00*sqrt(12.58),
                                   hatse.k=1/sqrt(12.58))  


    boundaries3 <- InterimAnalysis(SeqCR.object=boundaries2,
                                   hatbeta.k=-1.21*sqrt(21.11),
                                   hatse.k=1/sqrt(21.11))

    boundaries4 <- InterimAnalysis(SeqCR.object=boundaries3,
                                   hatbeta.k=-0.73*sqrt(30.55),
                                   hatse.k=1/sqrt(30.55))

    boundaries5 <- InterimAnalysis(SeqCR.object=boundaries4,
                                   hatbeta.k=-0.87*sqrt(33.28),
                                   hatse.k=1/sqrt(33.28))

    ## {{{ Maximal difference in boundaries
    Diffcomputed <- c(boundaries1$boundaries$Obs.b.k[1]-3.23,
                      boundaries2$boundaries$Obs.b.k[2]-2.76,
                      boundaries3$boundaries$Obs.b.k[3]-2.44,
                      boundaries4$boundaries$Obs.b.k[4]-2.16,
                      boundaries5$boundaries$Obs.a.k[5]-2.14)
    ## max(Diffcomputed)
    ## }}}       
    expect_equivalent(Diffcomputed,rep(0,length(Diffcomputed)),tolerance=5e-3)
  
})
## }}}


## {{{ Test reproduce values of Table 7.7 p 165 of Jennison and Turnbull's book (Expected sample size)


test_that("Expected sample size : Table 7.7 p. 165 : rho=2, K=4, beta=0.2, alpha=0.05",{
    MyResults <- planboundaries(rho=2,  
                                alpha=0.05,
                                beta=0.2,
                                K=4,         
                                theta=log(2),
                                abseps = 1e-08,
                                binding=TRUE,
                                direction="smaller",
                                Trace=FALSE,
                                nWhileMax=30,
                                toldiff=1e-06,
                                mycoefMax= 1.5,       
                                mycoefL=1             
                                )

    AverageN <- ExpectSampleSize(MyResults)

    expect_equivalent(AverageN$N[2],80.2,tolerance=5e-2)
})


test_that("Expected sample size : Table 7.7 p. 165 : rho=3, K=5, beta=0.2, alpha=0.05",{
    MyResults <- planboundaries(rho=3,  
                                alpha=0.05,
                                beta=0.2,
                                K=5,         
                                theta=log(2),
                                abseps = 1e-08,
                                binding=TRUE,
                                direction="smaller",
                                Trace=FALSE,
                                nWhileMax=30,
                                toldiff=1e-06,
                                mycoefMax= 1.5,       
                                mycoefL=1             
                                )

    AverageN <- ExpectSampleSize(MyResults)

    expect_equivalent(AverageN$N[3],77.8,tolerance=5e-2)
})


test_that("Expected sample size : Table 7.7 p. 165 : rho=3, K=1, beta=0.2, alpha=0.05",{
    MyResults <- planboundaries(rho=3,  
                                alpha=0.05,
                                beta=0.2,
                                K=1,         
                                theta=log(2),
                                abseps = 1e-08,
                                binding=TRUE,
                                direction="smaller",
                                Trace=FALSE,
                                nWhileMax=30,
                                toldiff=1e-06,
                                mycoefMax= 1.5,       
                                mycoefL=1             
                                )

    AverageN <- ExpectSampleSize(MyResults)

    expect_equivalent(AverageN$N[3],100,tolerance=5e-2)
})


test_that("Expected sample size : Table 7.8 p. 166 : rho=2, K=3, beta=0.1, alpha=0.05",{
    MyResults <- planboundaries(rho=2,  
                                alpha=0.05,
                                beta=0.1,
                                K=3,         
                                theta=log(1.5),
                                abseps = 1e-08,
                                binding=TRUE,
                                direction="smaller",
                                Trace=FALSE,
                                nWhileMax=30,
                                toldiff=1e-06,
                                mycoefMax= 1.5,       
                                mycoefL=1             
                                )

    AverageN <- ExpectSampleSize(MyResults)

    expect_equivalent(AverageN$N[2],83.9,tolerance=5e-2)
})

## }}}


## {{{ Test reproduce values of Fig 25.6 (Average sample size) of Chapter Jennison and Turnbull 2013 (and numbers in Section 25.5)

test_that("Reproduce values of Fig 25.6 of Chapter Jennison and Turnbull 2013 (and numbers in Section 25.5)",{
    MyResults <- planboundaries(rho=2,  
                                  alpha=0.025,
                                  beta=0.2,
                                  K=5,         
                                  theta=0.5,
                                  abseps = 1e-08,
                                  binding=TRUE,
                                  direction="smaller",
                                  Trace=FALSE,
                                  nWhileMax=30,
                                  toldiff=1e-06,
                                  mycoefMax= 1.2,       
                                  mycoefL=1             
                                  )

    AverageN <- ExpectSampleSize(MyResults, thetaValues = c(0,0.5,0.35/0.5,1))
    round(AverageN$N*(MyResults$If*4/100),1)

    # To reproduce the plot of Fig 25.6 of Chapter Jennison and Turnbull 2013
    ## AverageN <- ExpectSampleSize(MyResults,
    ## thetaValues = log(seq(from=exp(-0.5),to=exp(1),length.out=100))/MyResults$theta)
    ## plot(log(seq(from=exp(-0.5),to=exp(1),length.out=100)),
    ## type="l", lwd=2,
    ## AverageN$N*(MyResults$If*4/100),
    ## ylim=c(0,150),xlab="Log HR",ylab="E(number of events)")
    ## abline(h=(MyResults$If*4),lty=2, lwd=2,)
    ## legend("bottomright",c("fixed sample size (i.e. non sequential design)","sequential design"),lty=2:1,lwd=2,col="black")
    
    expect_equivalent(AverageN$N[c(1,3,4)]*(MyResults$If*4/100),c(72.9,100.9,94.5),tolerance=5e-2)
})


## }}}


## {{{ Comparison 1 with gsDesign

test_that("Comparison 1 with gsDesign",{

    ## {{{ Parameters
    Myrho <- 2
    Myalpha <- 0.025
    Mybeta <-0.2
    MyK <- 5
    Mytheta <- log(2)
    myvalues <- seq(from=0,to=1.5,by=0.1)
    mybinding <- TRUE
    ## }}}

    ## {{{ Estimates
    MyResults <- planboundaries(rho=Myrho,  
                                alpha=Myalpha,
                                beta=Mybeta,
                                K=MyK,         
                                theta=Mytheta,
                                abseps = 1e-07,
                                binding=mybinding,
                                direction="smaller",
                                Trace=TRUE,
                                nWhileMax=30,
                                toldiff=1e-06,
                                mycoefMax= 1.5,       
                                mycoefL=1             
                                )

    designCI <- gsDesign(k=MyK,
                         test.type = ifelse(mybinding,3,4), 
                         alpha=Myalpha, 
                         beta=Mybeta, 
                         delta=Mytheta,
                         sfu=sfPower,
                         sfl=sfPower,
                         sfupar=Myrho,
                         sflpar=Myrho)
    ## }}}

    ## {{{ Compare boundaries 
    bounds <- cbind(c(MyResults$boundaries$a.k,
                      MyResults$boundaries$b.k),
                    c(designCI$lower$bound,designCI$upper$bound))
    rownames(bounds) <- c(paste("lower, k==",1:MyResults$K),paste("upper, k=",1:MyResults$K))
    colnames(bounds) <- c("Mine","gsDesign")
    bounds <- cbind(bounds,abs(bounds[,"Mine"]-bounds[,"gsDesign"]))
    colnames(bounds)[3] <- "abs. Diff"
    bounds
    ## }}}

    ## {{{ Compare fixed sample size information
    bounds <- rbind(bounds,c(MyResults$If,designCI$n.fix,abs(MyResults$If-designCI$n.fix)))
    rownames(bounds)[nrow(bounds)] <- "If"
    bounds
    ## }}}

    ## {{{ Compare Average sample size

    avSSgsDesign <- gsProbability(d=designCI, theta=designCI$delta*myvalues)$en
    names(avSSgsDesign) <- paste("theta=",round(designCI$delta*myvalues,3),
                                 "(",myvalues*100,"%)",sep="")
    # This gives the expected number of events.
    # In the context of time to event setting with 1:1 randomization ratio,
    # we multiply by 4 to get the sample size (if there is no censoring).
    round(avSSgsDesign*4,1)


    ResAvNCont <- ExpectSampleSize(MyResults,
                                   actual.binding = TRUE, # always decide to continue even if we can stop for futility (if actual.binding = TRUE)
                                   thetaValues=myvalues
                                   ) 
    ResAvNCont$N[which(myvalues==1)] #  Expected sample size under H1: theta=delta (i.e. equal the value that we actually assume for computing the power)
    # In fact it is the expected information on termination, as in Jennison and Turnbull's book (see Par 7.3.3 p. 165),
    # as expressed in percentage of the corresponding fixed sample size.

    avSSmine <- ResAvNCont$N*(MyResults$If*4/100) 
    names(avSSmine) <- names(avSSgsDesign)

    bounds <- rbind(bounds,cbind(avSSmine,avSSgsDesign*4,abs(avSSmine-avSSgsDesign*4)))
    bounds
    ## }}}

    ## {{{ Maximal difference in boundaries and in expected sample size
    maxDevBoundaries <- max(bounds[1:(2*MyK+1),"abs. Diff"])*100
    maxDevExpSs <- max(bounds[(2*MyK+2):(2*MyK+1+length(myvalues)),"abs. Diff"])
    ## }}}    

   
    expect_equivalent(c(maxDevBoundaries,maxDevExpSs),c(0,0),tolerance=1e-2)
})


## }}}

## {{{ Comparison 2 with gsDesign

test_that("Comparison 2 with gsDesign",{

    ## {{{ Parameters
    Myrho <- 2
    Myalpha <- 0.05
    Mybeta <-0.1
    MyK <- 3
    Mytheta <- log(1.5)
    myvalues <- seq(from=0,to=1.75,by=0.1)
    mybinding <- FALSE    
    ## }}}

    ## {{{ Estimates
    MyResults <- planboundaries(rho=Myrho,  
                                alpha=Myalpha,
                                beta=Mybeta,
                                K=MyK,         
                                theta=Mytheta,
                                abseps = 1e-07,
                                binding=mybinding,
                                direction="smaller",
                                Trace=TRUE,
                                nWhileMax=30,
                                toldiff=1e-06,
                                mycoefMax= 1.5,       
                                mycoefL=1             
                                )

    designCI <- gsDesign(k=MyK,
                         test.type = ifelse(mybinding,3,4), 
                         alpha=Myalpha, 
                         beta=Mybeta, 
                         delta=Mytheta,
                         sfu=sfPower,
                         sfl=sfPower,
                         sfupar=Myrho,
                         sflpar=Myrho)
    ## }}}

    ## {{{ Compare boundaries 
    bounds <- cbind(c(MyResults$boundaries$a.k,
                      MyResults$boundaries$b.k),
                    c(designCI$lower$bound,designCI$upper$bound))
    rownames(bounds) <- c(paste("lower, k==",1:MyResults$K),paste("upper, k=",1:MyResults$K))
    colnames(bounds) <- c("Mine","gsDesign")
    bounds <- cbind(bounds,abs(bounds[,"Mine"]-bounds[,"gsDesign"]))
    colnames(bounds)[3] <- "abs. Diff"
    bounds
    ## }}}

    ## {{{ Compare fixed sample size information
    bounds <- rbind(bounds,c(MyResults$If,designCI$n.fix,abs(MyResults$If-designCI$n.fix)))
    rownames(bounds)[nrow(bounds)] <- "If"
    bounds
    ## }}}

    ## {{{ Compare Average sample size

    avSSgsDesign <- gsProbability(d=designCI, theta=designCI$delta*myvalues)$en
    names(avSSgsDesign) <- paste("theta=",round(designCI$delta*myvalues,3),
                                 "(",myvalues*100,"%)",sep="")
    # This gives the expected number of events.
    # In the context of time to event setting with 1:1 randomization ratio,
    # we multiply by 4 to get the sample size (if there is no censoring).
    round(avSSgsDesign*4,1)


    ResAvNCont <- ExpectSampleSize(MyResults,
                                   actual.binding = TRUE, # always decide to continue even if we can stop for futility (if actual.binding = TRUE)
                                   thetaValues=myvalues
                                   ) 
    ResAvNCont$N[which(myvalues==1)] #  Expected sample size under H1: theta=delta (i.e. equal the value that we actually assume for computing the power)
    # In fact it is the expected information on termination, as in Jennison and Turnbull's book (see Par 7.3.3 p. 165),
    # as expressed in percentage of the corresponding fixed sample size.

    avSSmine <- ResAvNCont$N*(MyResults$If*4/100) 
    names(avSSmine) <- names(avSSgsDesign)

    bounds <- rbind(bounds,cbind(avSSmine,avSSgsDesign*4,abs(avSSmine-avSSgsDesign*4)))
    bounds
    ## }}}

    ## {{{ Maximal difference in boundaries and in expected sample size
    maxDevBoundaries <- max(bounds[1:(2*MyK+1),"abs. Diff"])*100
    maxDevExpSs <- max(bounds[(2*MyK+2):(2*MyK+1+length(myvalues)),"abs. Diff"])
    ## }}}    

   
    expect_equivalent(c(maxDevBoundaries,maxDevExpSs),c(0,0),tolerance=1e-2)
})


## }}}

## {{{ Comparison 3 with gsDesign (Neotrans)

test_that("Comparison 3 with gsDesign (Neotrans))",{

    ## {{{ parameters
    myrho <- 1:3
    myalpha <- 0.05
    mybeta <- 0.2
    ## myK <- 4
    myK <- 2
    mytheta <- log(2)
    mybinding <- FALSE
    mydirection <- "smaller"
    myvalues <- seq(from=0,to=1.5,by=0.05)
    ## }}}


    
    ## {{{  My estimates
    MyResults1 <- planboundaries(rho=myrho[1],  
                                 alpha=myalpha,
                                 beta=mybeta,
                                 K=myK,         
                                 theta=mytheta,
                                 abseps = 1e-06,
                                 binding=mybinding,
                                 direction=mydirection,
                                 Trace=FALSE,
                                 nWhileMax=30,
                                 toldiff=1e-05,
                                 mycoefMax= 1.5,       
                                 mycoefL=1             
                                 )

    MyResults2 <- planboundaries(rho=myrho[2],  
                                 alpha=myalpha,
                                 beta=mybeta,
                                 K=myK,         
                                 theta=mytheta,
                                 abseps = 1e-08,
                                 binding=mybinding,
                                 direction=mydirection,
                                 Trace=FALSE,
                                 nWhileMax=30,
                                 toldiff=1e-06,
                                 mycoefMax= 1.5,       
                                 mycoefL=1             
                                 )
    MyResults3 <- planboundaries(rho=myrho[3],  
                                 alpha=myalpha,
                                 beta=mybeta,
                                 K=myK,         
                                 theta=mytheta,
                                 abseps = 1e-08,
                                 binding=mybinding,
                                 direction=mydirection,
                                 Trace=FALSE,
                                 nWhileMax=30,
                                 toldiff=1e-06,
                                 mycoefMax= 1.5,       
                                 mycoefL=1             
                                 )
    
    ResAvNCont1 <- ExpectSampleSize(MyResults1,
                                    actual.binding = TRUE,
                                    thetaValues=myvalues
                                    ) 
                                    
    ResAvNCont2 <- ExpectSampleSize(MyResults2,
                                    actual.binding = TRUE,
                                    thetaValues=myvalues
                                    ) 

    ResAvNCont3 <- ExpectSampleSize(MyResults3,
                                    actual.binding = TRUE,
                                    thetaValues=myvalues
                                    )

    allSs <- cbind(ResAvNCont1$N*(MyResults1$If*4)/100,
                   ResAvNCont2$N*(MyResults2$If*4)/100,
                   ResAvNCont3$N*(MyResults3$If*4)/100)
    colnames(allSs) <- paste("rho=",1:3,sep="")
    round(allSs,1)
    ## }}}


    ## {{{ gsDesign
    designCI1 <- gsDesign(k=myK,
                          test.type = ifelse(mybinding,3,4), 
                          alpha=myalpha, 
                          beta=mybeta, 
                          delta=mytheta,
                          sfu=sfPower,
                          sfl=sfPower,
                          sfupar=myrho[1],
                          sflpar=myrho[1])
    designCI2 <- gsDesign(k=myK,
                          test.type = ifelse(mybinding,3,4), 
                          alpha=myalpha, 
                          beta=mybeta, 
                          delta=mytheta,
                          sfu=sfPower,
                          sfl=sfPower,
                          sfupar=myrho[2],
                          sflpar=myrho[2])
    designCI3 <- gsDesign(k=myK,
                          test.type = ifelse(mybinding,3,4), 
                          alpha=myalpha, 
                          beta=mybeta, 
                          delta=mytheta,
                          sfu=sfPower,
                          sfl=sfPower,
                          sfupar=myrho[3],
                          sflpar=myrho[3])

    avSSgsDesign1 <- gsProbability(d=designCI1, theta=designCI1$delta*myvalues)$en
    avSSgsDesign2 <- gsProbability(d=designCI2, theta=designCI2$delta*myvalues)$en
    avSSgsDesign3 <- gsProbability(d=designCI3, theta=designCI3$delta*myvalues)$en
    
    allSsgsDesign <- cbind(avSSgsDesign1,
                           avSSgsDesign2,
                           avSSgsDesign3)*4
    rownames(allSsgsDesign) <- myvalues

    colnames(allSsgsDesign) <- paste("rho=",1:3,sep="")
        
    round(allSsgsDesign,1)
    ## }}}

    ## {{{ Comparison
    Tocompare <- cbind(allSs,allSsgsDesign,abs(allSs-allSsgsDesign))
    colnames(Tocompare) <- c(paste("Mine rho=",1:3,sep=""),
                             paste("gsDesing rho=",1:3,sep=""),
                             paste("abs diff rho=",1:3,sep=""))
    Tocompare
    ## }}}
    
    expect_equivalent(apply(Tocompare[,7:9],2,max),rep(0,3),tolerance=5e-3)
})
    
## }}}
