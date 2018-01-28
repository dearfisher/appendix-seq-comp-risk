## {{{ Description
# This function computes expected sampe sizes in the case of group-sequential clinical trial (which are planned with the function planboundaries)

# In fact it computes the expected information on termination, as in
# Jennison and Turnbull's book (see Par 7.3.3 p. 165), as expressed in
# percentage of the corresponding fixed sample size.
# It is then easy to convert into average sample sizes (when ignoring some minor roundings)
# by multiplying by the fixed sample size (or e.g. the number of events with time to event data) and
# dividing by 100. this is because the information is proportional to the sample size (or to the number of events).
# The computation assumes that the levels of information observed at each analysis
# match with those which have been planned.


## Input
# SeqCR.object: an object fitted from function planboundaries
# thetaValues = c(0,0.5,1) fraction of theta values (i.e. expected differnces) for the effect size under H1 (default are c(0,0.5,1), as in Jennison and Turnbull's book Table 7.7 p 165)
# abseps = 1e-06 : tolerance for precision when computing integrals
# actual.binding = FALSE   # if SeqCR.object$binding=FALSE, we plan boundaries such that we can, but we do not have  to, continue. In this case, choosing actual.binding = FALSE  will compute the average sample size as if we always decide to continue even if we can stop for futility. Choosing actual.binding = TRUE  will compute the average sample size as we decide to stop for futility as early as we can.
## }}}

## {{{ Function

ExpectSampleSize <- function(SeqCR.object,
                             thetaValues = c(0,0.5,1),
                             actual.binding = FALSE,  
                             abseps = 1e-06){
    ## {{{ preliminaries
    x <- SeqCR.object$boundaries
    K <- SeqCR.object$K
    thea <- x$a.k
    theb <- x$b.k
    direction <- SeqCR.object$direction
    Imax <- SeqCR.object$Imax
    theta <- SeqCR.object$theta

    binding <- SeqCR.object$binding    
    if(!binding){
        binding <- actual.binding
    }
    ## print(paste("binding is ",binding))
    ## }}}
    ## {{{ compute variance-covariance matrix of vector (Z_1,...,Z_k)
    # assuming equally spaced information
    sigmaZk <- diag(1,K)
    for(i in 1:K){
        for(j in i:K){
            sigmaZk[i,j] <- sqrt(i/j)
            sigmaZk[j,i] <- sqrt(i/j)
        }
    }
    ## }}}
    ## {{{ initialize probabilities of stopping at analysis k    
    ProbaStop <- matrix(NA,nrow=K,ncol=length(thetaValues))
    colnames(ProbaStop) <- thetaValues
    rownames(ProbaStop) <- paste("k=",1:K,sep="")
    ## }}}
    ## {{{ compute the mean of the multivariate normal distribution under the alternative H1
    # assuming equally spaced information
    ## browser()
    thetheta <- matrix(theta*sqrt(((1:K)/K)*Imax),ncol=length(thetaValues),nrow=K)*matrix(thetaValues,ncol=length(thetaValues),nrow=K,byrow=TRUE)
    colnames(thetheta) <- thetaValues
    rownames(thetheta) <- paste("k=",1:K,sep="")
    ## }}}
    ## {{{  compute probabilities of stopping at analysis k
    ## {{{ k=1
    if(direction=="smaller"){
        for(j in 1:length(thetheta[1,])){
            if(binding){
                ## browser()
                ProbaStop[1,j] <- 1-pnorm(q=theb[1],sd=1,mean=thetheta[1,j]) + pnorm(q=thea[1],sd=1,mean=thetheta[1,j]) # Proba of stopping for efficacy + Proba of stopping for futility
            }else{
                ProbaStop[1,j] <- 1-pnorm(q=theb[1],sd=1,mean=thetheta[1,j]) # Proba of stopping for efficacy 
            }
        }
    }else{
        for(j in 1:length(thetheta[1,])){
            if(binding){            
                ProbaStop[1,j] <- pnorm(q=theb[1],sd=1,mean=thetheta[1,j]) + 1-pnorm(q=thea[1],sd=1,mean=thetheta[1,j]) # Proba of stopping for efficacy or futility
            }else{
                ProbaStop[1,j] <- pnorm(q=theb[1],sd=1,mean=thetheta[1,j])  # Proba of stopping for efficacy               
            }
        }
    }
    ## print("1 : DONE")
    ## }}}
    ## {{{ k>1
    if(K>1){
        if(direction=="smaller"){
            ## print("direction==smaller")
            ## {{{ 
            for(k in 2:K){
                for(j in 1:length(thetheta[k,])){
                    ## print(paste("j=",j))
                    if(binding){
                        ## print("binding==TRUE")
                        if(k!=K){
                            ## print("k!=K")                            
                            ## {{{ 
                            # Proba of stopping for efficacy + Proba of stopping for futility
                            ProbaStop[k,j] <- pmvnorm(lower = c(thea[1:(k-1)],theb[k]),  
                                                      upper = c(theb[1:(k-1)],Inf),
                                                      mean = thetheta[1:k,j],
                                                      sigma= sigmaZk[1:k,1:k],
                                                      abseps = abseps)  
                            ProbaStop[k,j] <- ProbaStop[k,j] + pmvnorm(lower = c(thea[1:(k-1)],-Inf),  
                                                                       upper = c(theb[1:(k-1)],thea[k]),
                                                                       mean = thetheta[1:k,j],
                                                                       sigma= sigmaZk[1:k,1:k],
                                                                       abseps = abseps)
                            ## }}}
                        }else{                            
                            ## print("k==K")                            
                            ## {{{ 
                            ProbaStop[k,j] <- pmvnorm(lower = thea[1:(k-1)],  
                                                      upper = theb[1:(k-1)],
                                                      mean = thetheta[1:(k-1),j],
                                                      sigma= sigmaZk[1:(k-1),1:(k-1)],
                                                      abseps = abseps)
                            ## print(ProbaStop[k,j])
                            ## }}}
                        }
                    }else{
                        ## print("binding==FALSE")
                        if(k!=K){
                            ## print("k!=K")                            
                            ## {{{ 
                            # Proba of stopping for efficacy 
                            ProbaStop[k,j] <- pmvnorm(lower = c(rep(-Inf,k-1),theb[k]),
                                                      upper = c(theb[1:(k-1)],Inf),
                                                      mean = thetheta[1:k,j],
                                                      sigma= sigmaZk[1:k,1:k],
                                                      abseps = abseps)
                            ## print(ProbaStop[k,j])
                            ## }}}
                        }else{
                            ## {{{
                            ## print("k==K")                            
                            ProbaStop[k,j] <- pmvnorm(lower = rep(-Inf,k-1),  
                                                      upper = theb[1:(k-1)],
                                                      mean = thetheta[1:(k-1),j],
                                                      sigma= sigmaZk[1:(k-1),1:(k-1)],
                                                      abseps = abseps)
                            ## print(ProbaStop[k,j])
                            ## }}}
                        }
                    }
                }
            }
            ## }}}
        }else{
            ## {{{ 
            for(k in 2:K){
                for(j in 1:length(thetheta[k,])){
                    if(binding){
                        if(k!=K){
                            ## {{{ 
                            # Proba of stopping for efficacy + Proba of stopping for futility
                            ProbaStop[k,j] <- pmvnorm(lower = c(theb[1:(k-1)],-Inf),
                                                      upper = c(thea[1:(k-1)],theb[k]),
                                                      mean = thetheta[1:k,j],
                                                      sigma= sigmaZk[1:k,1:k],
                                                      abseps = abseps)
                            ProbaStop[k,j] <- ProbaStop[k,j] + pmvnorm(lower = c(theb[1:(k-1)],thea[k]),
                                                                       upper = c(thea[1:(k-1)],Inf),
                                                                       mean = thetheta[1:k,j],
                                                                       sigma= sigmaZk[1:k,1:k],
                                                                       abseps = abseps)
                            ## }}}
                        }else{
                            ## {{{ 
                            ProbaStop[k,j] <- pmvnorm(lower = theb[1:(k-1)],
                                                      upper = thea[1:(k-1)],
                                                      mean = thetheta[1:k,j],
                                                      sigma= sigmaZk[1:k,1:k],
                                                      abseps = abseps)
                            ## }}}
                        }
                    }else{
                        if(k!=K){
                            ## {{{                            
                            # Proba of stopping for efficacy 
                            ProbaStop[k,j] <- pmvnorm(lower = c(theb[1:(k-1)],-Inf),
                                                      upper = c(rep(Inf,k-1),theb[k]),
                                                      mean = thetheta[1:k,j],
                                                      sigma= sigmaZk[1:k,1:k],
                                                      abseps = abseps)
                            ## }}}
                        }else{
                            ## {{{ 
                            ProbaStop[k,j] <- pmvnorm(lower = theb[1:(k-1)],
                                                      upper = rep(Inf,k-1),
                                                      mean = thetheta[1:k,j],
                                                      sigma= sigmaZk[1:k,1:k],
                                                      abseps = abseps)
                            ## }}}
                        }    
                    }
                }
            }
            ## }}}
        }
    }
    ## }}}
    ## }}}
    ## {{{ Create output
    N <- colSums(ProbaStop*matrix((1:K)/K,
                                  ncol=length(thetaValues),
                                  nrow=K,byrow=FALSE)
                 )*SeqCR.object$coef*100
    out <- list(N=N,Prob=ProbaStop)
    ## }}}
    out
}

## }}}
