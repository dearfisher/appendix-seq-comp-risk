InterimAnalysis <- function(SeqCR.object,   # object of class "SeqCR", i.e. as returned by planboundaries or InterimAnalysis
                            k=NULL,         # for which interim analysis we want to add the results
                            hatbeta.k,      # estimated effect (on the scale of the test statistc for which If and Imax relate to one over the variance, e.g. estimated log(HR))
                            hatse.k         # estimated s.e. of the effect
                            ){
    ## {{{ Add column if needed
    x <- SeqCR.object
    d <- x$boundaries
    ##
    if(!("Obs.I.k" %in% names(d)) ){
        d$Obs.I.k <- rep(NA,nrow(d))
    }
    if(!("Obs.z.k" %in% names(d)) ){
        d$Obs.z.k <- rep(NA,nrow(d))
    }
    if(!("Obs.a.k" %in% names(d)) ){
        d$Obs.a.k <- rep(NA,nrow(d))
    }
    if(!("Obs.b.k" %in% names(d)) ){
        d$Obs.b.k <- rep(NA,nrow(d))
    }
    if(!("Obs.Inc.Type.II.error" %in% names(d)) ){
        d$Obs.Inc.Type.II.error <- rep(NA,nrow(d))
    }
    if(!("Obs.Inc.Type.I.error" %in% names(d)) ){
        d$Obs.Inc.Type.I.error <- rep(NA,nrow(d))
    }
    if(!("Obs.Type.I.error" %in% names(d)) ){
        d$Obs.Type.I.error <- rep(NA,nrow(d))
    }
    if(!("Obs.Type.II.error" %in% names(d)) ){
        d$Obs.Type.II.error <- rep(NA,nrow(d))
    }
    ## }}}
    ## {{{ if k is NULL 
    if(is.null(k)){
        k <- min(which(is.na(d$Obs.z.k)))
    }
    ## }}}

    ## {{{ add results of interim analysis
    ## browser()
    d$Obs.z.k[k] <- hatbeta.k/hatse.k
    d$Obs.I.k[k] <- 1/(hatse.k)^2
    d$Obs.Type.I.error[k] <- f(I=d$Obs.I.k[k],rho=x$rho,alpha=x$alpha,Imax=x$Imax)
    d$Obs.Type.II.error[k] <- g(I=d$Obs.I.k[k],rho=x$rho,beta=x$beta,Imax=x$Imax)   
    if(k==1){
        d$Obs.Inc.Type.I.error[1] <- d$Obs.Type.I.error[1]    # - 0
        d$Obs.Inc.Type.II.error[1] <- d$Obs.Type.II.error[1]  # - 0
    }
    if(k>1){
        d$Obs.Inc.Type.I.error[k] <- d$Obs.Type.I.error[k] - d$Obs.Type.I.error[k-1]
        d$Obs.Inc.Type.II.error[k] <- d$Obs.Type.II.error[k] - d$Obs.Type.II.error[k-1]
    }
    ## }}}    
    ## {{{ Compute new boundaries
    thetheta <- SeqCR.object$theta*sqrt(d$Obs.I.k)
    ## {{{ k=1
    # for k=1 it does not matter whether binding is TRUE or FALSE
    if(k==1){
        if(SeqCR.object$direction=="smaller"){
            d$Obs.b.k[1] <- qnorm(p=1-d$Obs.Inc.Type.I.error[1],mean=0,sd=1)         # compute under the null (Ho)
            d$Obs.a.k[1] <- qnorm(p=d$Obs.Inc.Type.II.error[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)
            # handle over running
            if(d$Obs.b.k[1]<d$Obs.a.k[1]){
                d$Obs.a.k[1] <- d$Obs.b.k[1]
                d$Obs.Inc.Type.II.error[1] <- pnorm(d$Obs.a.k[1],mean=thetheta[1],sd=1)
            }
        }else{
            ## i.e if SeqCR.object$direction=="greater"
            d$Obs.b.k[1] <- qnorm(p=d$Obs.Inc.Type.I.error[1],mean=0,sd=1)             # compute under the null (Ho)
            d$Obs.a.k[1] <- qnorm(p=1-d$Obs.Inc.Type.II.error[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)
            # handle over running
            if(d$Obs.b.k[1] > d$Obs.a.k[1]){
                d$Obs.a.k[1] <- d$Obs.b.k[1]
                d$Obs.Inc.Type.II.error[1] <- 1-pnorm(d$Obs.a.k[1],mean=thetheta[1],sd=1)
            }
        }
    }

    ## }}}
    ## {{{ k>1
    if(k>1){
        ## {{{ compute variance-covariance matrix
        # of vector (Z_1,...,Z_k) assuming equally spaced information
        sigmaZk <- diag(1,k)
        for(i in 1:k){
            for(j in i:k){
                sigmaZk[i,j] <- sqrt(d$Obs.I.k[i]/d$Obs.I.k[j])
                sigmaZk[j,i] <- sqrt(d$Obs.I.k[i]/d$Obs.I.k[j])
            }
        }
        ## }}}
        if(SeqCR.object$binding){
            ## {{{ 
            if(SeqCR.object$direction=="smaller"){
                ## {{{ b_k by solving what follows 
                ## ## theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows (when binding = TRUE )
                ## try(
                d$Obs.b.k[k] <- uniroot(function(x){pmvnorm(lower = c(d$Obs.a.k[1:(k-1)],x),
                                                            upper = c(d$Obs.b.k[1:(k-1)],Inf),
                                                            mean=rep(0,k),
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.I.error[k]},
                                        lower = d$Obs.a.k[1],
                                        upper = d$Obs.b.k[1],
                                        tol = SeqCR.object$abseps)$root
                ## , silent = TRUE)
                ## IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
                ## }}}
                ## {{{ a_k by solving what follows 
                ## if(IsbkOK){
                ## thea[k] <- theb[k] # just to handle cases in which there is no root in what follows 
                ## try(
                d$Obs.a.k[k] <- uniroot(function(x){pmvnorm(lower = c(d$Obs.a.k[1:(k-1)],-Inf),
                                                            upper = c(d$Obs.b.k[1:(k-1)],x),
                                                            mean=thetheta[1:k],
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.II.error[k]},
                                        lower = d$Obs.a.k[1], 
                                        upper = d$Obs.b.k[1], 
                                        tol = SeqCR.object$abseps)$root
                ## , silent = TRUE)
                ## }else{
                ## thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
                ## }
                ## }}}
                ## {{{ handle either over running or under running
                if(d$Obs.b.k[k]<d$Obs.a.k[k] | k==SeqCR.object$K){
                    d$Obs.a.k[k] <- d$Obs.b.k[k]
                    d$Obs.Inc.Type.II.error[k] <- pmvnorm(lower = c(d$Obs.a.k[1:(k-1)],-Inf),
                                                          upper = c(d$Obs.b.k[1:(k-1)],d$Obs.a.k[k]),
                                                          mean=thetheta[1:k],
                                                          sigma= sigmaZk,
                                                          abseps = SeqCR.object$abseps)
                    d$Obs.Type.II.error[k] <- d$Obs.Type.II.error[k-1] + d$Obs.Inc.Type.II.error[k]
                }
                ## }}}
            }
            if(SeqCR.object$direction=="greater"){
                ## {{{ b_k by solving what follows 
                ## theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows  (when binding = TRUE )
                ## try(
                d$Obs.b.k[k] <- uniroot(function(x){pmvnorm(lower = c(d$Obs.b.k[1:(k-1)],-Inf),
                                                            upper = c(d$Obs.a.k[1:(k-1)],x),
                                                            mean=rep(0,k), 
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.I.error[k]},
                                        lower = d$Obs.b.k[1],
                                        upper = d$Obs.a.k[1],
                                        tol = SeqCR.object$abseps)$root
                ## , silent = TRUE)
                ## IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
                ## }}}
                ## {{{ a_k by solving what follows 
                ## if(IsbkOK){
                ## thea[k] <- theb[k] # just to handle cases in which there is no root in what follows                            
                ## try(
                d$Obs.a.k[k] <- uniroot(function(x){pmvnorm(lower = c(d$Obs.b.k[1:(k-1)],x),
                                                            upper = c(d$Obs.a.k[1:(k-1)],Inf),
                                                            mean=thetheta[1:k],
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.II.error[k]},
                                        lower = d$Obs.b.k[1],
                                        upper = d$Obs.a.k[1],
                                        tol = SeqCR.object$abseps)$root
                ## , silent = TRUE)
                ## }else{
                ## thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
                ## }
                ## }}}
                ## {{{ handle either over running or under running
                if(d$Obs.b.k[k]>d$Obs.a.k[k] | k==SeqCR.object$K){
                    d$Obs.a.k[k] <- d$Obs.b.k[k]
                    d$Obs.Inc.Type.II.error[k] <- pmvnorm(lower = c(d$Obs.b.k[1:(k-1)],d$Obs.a.k[k]),
                                                          upper = c(d$Obs.a.k[1:(k-1)],Inf),
                                                          mean=thetheta[1:k],
                                                          sigma= sigmaZk,
                                                          abseps = SeqCR.object$abseps)
                    d$Obs.Type.II.error[k] <- d$Obs.Type.II.error[k-1] + d$Obs.Inc.Type.II.error[k]
                }
                ## }}}
            }
            ## }}}
        }else{
            ## {{{ 
            if(SeqCR.object$direction=="smaller"){
                ## {{{ b_k by solving what follows 
                d$Obs.b.k[k] <- uniroot(function(x){pmvnorm(lower = c(rep(-Inf, k-1),x),
                                                            upper = c(d$Obs.b.k[1:(k-1)],Inf),
                                                            mean=rep(0,k),
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.I.error[k]},
                                        lower = d$Obs.a.k[1],
                                        upper = d$Obs.b.k[1],
                                        tol = SeqCR.object$abseps)$root
                ## }}}
                ## {{{ a_k by solving what follows 
                d$Obs.a.k[k] <- uniroot(function(x){pmvnorm(lower = c(d$Obs.a.k[1:(k-1)],-Inf),
                                                            upper = c(d$Obs.b.k[1:(k-1)],x),
                                                            mean=thetheta[1:k],
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.II.error[k]},
                                        lower = d$Obs.a.k[1], 
                                        upper = d$Obs.b.k[1], 
                                        tol = SeqCR.object$abseps)$root
                ## }}}
                ## {{{ handle either over running or under running
                if(d$Obs.b.k[k]<d$Obs.a.k[k] | k==SeqCR.object$K){
                    d$Obs.a.k[k] <- d$Obs.b.k[k]
                    d$Obs.Inc.Type.II.error[k] <- pmvnorm(lower = c(d$Obs.a.k[1:(k-1)],-Inf),
                                                          upper = c(d$Obs.b.k[1:(k-1)],d$Obs.a.k[k]),
                                                          mean=thetheta[1:k],
                                                          sigma= sigmaZk,
                                                          abseps = SeqCR.object$abseps)
                    d$Obs.Type.II.error[k] <- d$Obs.Type.II.error[k-1] + d$Obs.Inc.Type.II.error[k]
                }
                ## }}}
            }
            if(SeqCR.object$direction=="greater"){
                ## {{{ b_k by solving what follows 
                d$Obs.b.k[k] <- uniroot(function(x){pmvnorm(lower = d$Obs.b.k[1:k],
                                                            upper = c(rep(-Inf, k-1),x),
                                                            mean=rep(0,k), 
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.I.error[k]},
                                        lower = d$Obs.b.k[1],
                                        upper = d$Obs.a.k[1],
                                        tol = SeqCR.object$abseps)$root
                ## }}}
                ## {{{ a_k by solving what follows 
                d$Obs.a.k[k] <- uniroot(function(x){pmvnorm(lower = c(d$Obs.b.k[1:(k-1)],x),
                                                            upper = c(d$Obs.a.k[1:(k-1)],Inf),
                                                            mean=thetheta[1:k],
                                                            sigma= sigmaZk,
                                                            abseps = SeqCR.object$abseps) - d$Obs.Inc.Type.II.error[k]},
                                        lower = d$Obs.b.k[1],
                                        upper = d$Obs.a.k[1],
                                        tol = SeqCR.object$abseps)$root
                ## }}}
                ## {{{ handle either over running or under running
                if(d$Obs.b.k[k]<d$Obs.a.k[k] | k==SeqCR.object$K){
                    d$Obs.a.k[k] <- d$Obs.b.k[k]
                    d$Obs.Inc.Type.II.error[k] <- pmvnorm(lower = c(d$Obs.b.k[1:(k-1)],d$Obs.a.k[k]),
                                                          upper = c(d$Obs.a.k[1:(k-1)],Inf),
                                                          mean=thetheta[1:k],
                                                          sigma= sigmaZk,
                                                          abseps = SeqCR.object$abseps)
                    d$Obs.Type.II.error[k] <- d$Obs.Type.II.error[k-1] + d$Obs.Inc.Type.II.error[k]
                }
                ## }}}
            }
            ## }}}
        }
    }
    ## }}}        
    ## }}}
    ## {{{ Reorder columns
    d <- d[,c("a.k",
              "b.k",
              "Type.I.Error",
              "Type.II.Error",
              "Inc.Type.I",
              "Inc.Type.II",
              "I.k",                   
              "Obs.a.k",
              "Obs.b.k",
              "Obs.Type.I.error",
              "Obs.Type.II.error",
              "Obs.Inc.Type.I.error",
              "Obs.Inc.Type.II.error",                     
              "Obs.I.k",
              "Obs.z.k")]
    
    ## }}}
    ## {{{ create output
    x$boundaries <- d
    x
    ## }}}
}
