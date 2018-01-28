planboundaries <- function(rho=1,                # rho parameter of the rho-family spending functions (Kim-DeMets) for alpha and beta 
                           alpha=0.05,           # Type-I error (overall)
                           beta=0.2,             # Type-II error (overall)
                           K,                    # number of planned analyses
                           Imax=NULL,            # Imax, i.e. maximum information needed for given beta (type II-error), theta (expected difference), alpha (type I-error) and K. It can be given it is known. Otherwise it is computed from the  values given for alpha, beta, theta and K.
                           theta=0,              # expected effect under the alternative (should be on the scale of the test statistc for which If and Imax relate to one over the variance, e.g. theta=expected log(Hazard ratio))
                           abseps = 1e-06,       # tolerance for precision when finding roots or computing integrals
                           binding=TRUE,         # binding or non-binding futility boundary (i.e FALSE if it is allowed to continue after crossing the futility boundary)
                           direction="smaller",  # greater is for Ho= theta > 0, "smaller" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider smaller)
                           Trace=FALSE,          # Used only if Imax=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Imax (from  alpha, beta, theta and K).
                           nWhileMax=30,         # Used only if Imax=NULL. Maximum number of steps in the (root finding) algorithm to compute Imax (from  alpha, beta, theta and K)
                           toldiff= 1e-05,       # Used only if Imax=NULL. Maximum tolerated difference between lower and upper bounds at anaylis K (which souhld be zero), in the root finding algorithm, to find the value of Imax
                           tolcoef= 1e-04,       # Used only if Imax=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Imax=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if K=1)
                           mycoefMax= 1.2,       # Used only if Imax=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Imax=coef*If (in the root finding algorithm).
                           mycoefL=1,            # Used only if Imax=NULL. Lower limit of the interval (see mycoefMax)
                           myseed=2902           # seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
                           ){
    ## browser()
    ## {{{ set seed
    old <- .Random.seed # to save the current seed
    set.seed(myseed)
    ## }}}
    ## {{{ preliminaries
    mycoef <- NULL # initialize as needed for output
    if(direction=="smaller"){
        thea <- rep(-Inf,K)
        theb <- rep(Inf,K)     
    }else{
        if(direction=="greater"){
            thea <- rep(Inf,K)
            theb <- rep(-Inf,K)
        }else{
            stop("direction should be either greater or smaller")}
    }
    ## browser()
    if( (direction=="smaller" & theta<0) | (direction=="greater" & theta>0)){
        stop("The values given for arguments direction and theta are inconsistent.\n When direction=smaller, theta should be positive.\n When direction=greater, theta should be negative.")
    }
    # initialize
    thealpha <- rep(0,K)  # alpha spent up to step k
    thebeta <- rep(0,K)   # beta spent up to step k
    IncAlpha <- rep(0,K)  # alpha spent at step k
    IncBeta <- rep(0,K)   # beta spent at step k     
    # compute variance-covariance matrix of vector (Z_1,...,Z_k) assuming equally spaced information
    sigmaZk <- diag(1,K)
    for(i in 1:K){
        for(j in i:K){
            sigmaZk[i,j] <- sqrt(i/j)
            sigmaZk[j,i] <- sqrt(i/j)
        }
    }
    # compute If (see Jennison book page 87
    If <- (qnorm(1-alpha)+qnorm(1-beta))^2/theta^2
    if(Trace){
        cat("\n If computed as =",If,"\n")
    }
    ## }}}
    if(is.null(Imax)){
        ## {{{ Compute Imax from If and the other arguments (Recursive calls to the function)
        if(Trace){
            cat("\n We start the search of the value for coef=Imax/If. \n \n")
        }        
        ## {{{ initialize key values to be updated in the following loop
        nwhile <- 0
        mycoefL0 <- mycoefL
        mycoefU <- mycoefMax
        mycoef <- mycoefU
        ## }}}
        ## {{{ Is the interval within which to search for coef large enough ?      
        if(Trace){
            cat("\n Check whether the interval within which to search for coef large enough. \n")
        }
        xx <- planboundaries(rho=rho,
                             alpha=alpha,
                             beta=beta,
                             K=K,         
                             Imax=If*mycoefU,
                             theta=theta,
                             abseps=abseps,
                             toldiff=toldiff,
                             binding=binding,
                             direction=direction,
                             Trace=FALSE)
        thediff <- abs(xx$boundaries[K,"b.k"]-xx$boundaries[K,"a.k"])
        ## }}}       
        if(thediff==0){
            ## {{{ if yes, we search coef
            mycoef <- (mycoefL + mycoefU)/2
            thediff <- 2*toldiff
            if(Trace){
                cat("\n we start the search within [",mycoefL,",",mycoefU,"] \n")
            }
            while(nwhile < nWhileMax & thediff>toldiff & abs(mycoefL-mycoefU)> tolcoef){
                nwhile <- nwhile + 1
                if(Trace){
                    cat("\n Step :",nwhile,"(out of max.", nWhileMax,")")
                }
                xx <- planboundaries(rho=rho,
                                     alpha=alpha,
                                     beta=beta,
                                     K=K,         
                                     Imax=If*mycoef,
                                     theta=theta,
                                     abseps=abseps,
                                     toldiff=toldiff,
                                     binding=binding,
                                     direction=direction)
                thediff <- abs(xx$boundaries[K,"b.k"]-xx$boundaries[K,"a.k"])
                if(thediff>toldiff){
                    if(Trace){
                        cat("\n Value coef=",mycoef,"is too small  \n")
                        cat("\n coef=",mycoef,"leads to b.K-a.K=",thediff, "(whereas tol=",toldiff,") \n")
                        ## cat("\n b.K=",xx$boundaries[K,"b.k"]," and a.K=",xx$boundaries[K,"a.k"]," \n")
                    }
                    mycoefL <- (mycoefL+mycoefU)/2
                    if(Trace){
                        cat("\n we update the interval : [",mycoefL,",",mycoefU,"] \n")
                    }
                    mycoef <- (mycoefL+mycoefU)/2
                }
                if(thediff==0){
                    if(Trace){
                        cat("\n Value coef=",mycoef,"is too large  \n")
                        cat("\n coef=",mycoef,"leads to b.K-a.K=",thediff, "(whereas tol=",toldiff,") \n")
                        ## cat("\n b.K=",xx$boundaries[K,"b.k"]," and a.K=",xx$boundaries[K,"a.k"]," \n")
                    }
                    mycoefU <- (mycoefL+mycoefU)/2
                    if(Trace){
                        cat("\n we update the interval : [",mycoefL,",",mycoefU,"] \n")
                    }
                    mycoef <- (mycoefL+mycoefU)/2
                    thediff <- 2*toldiff
                }
                if((thediff<=toldiff & thediff!=0) | abs(mycoefL-mycoefU)<= tolcoef ){
                    if(Trace){
                        cat("\n coef value FOUND : coef=",mycoef,"\n (leads to b.K-a.K=",thediff, " and tol.=",toldiff," and search interval length is=",abs(mycoefL-mycoefU),"and tol.=",tolcoef,")\n")
                    }
                    Imax <- mycoef*If
                    if(Trace){
                        cat("\n Imax computed as=",Imax,"\n")
                    }
                    ## browser()
                    ## print("Imax created")
                }else{
                    if(nwhile==nWhileMax){
                        stop("Imax could not be computed presicely enough : we need to allow for more iterations in the algorithm : you should probably call the function again with a larger value for nWhileMax.")
                    }
                }
            }
            ## }}}
        }else{
            ## {{{ if no, we stop and explain why
            stop("The interval [mycoefL,mycoefMax]= [",mycoefL0,",",mycoefMax,"] is too small. You should probably call the function again with a larger value for mycoefMax and/or a lower (value >=1) for mycoefL \n")
            ## }}}
        }        

        ## }}}
    }else{
        mycoef <- Imax/If
    }
    ## browser()
    ## {{{ compute vector of means under the alternative H1
    # compute the mean of the multivariate normal distribution under the alternative H1 assuming equally spaced information
    ## browser()
    thetheta <- theta*sqrt(((1:K)/K)*Imax)
    ## }}} 
    ## {{{ case k=1 
    IncAlpha[1] <- f(I=(1/K)*Imax,rho=rho,alpha=alpha,Imax=Imax)
    IncBeta[1] <-  g(I=(1/K)*Imax,rho=rho,beta=beta,Imax=Imax)
    if(direction=="smaller"){
        theb[1] <- qnorm(p=1-IncAlpha[1],mean=0,sd=1)         # compute under the null (Ho)
        thea[1] <- qnorm(p=IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)
    }else{
        if(direction=="greater"){
            theb[1] <- qnorm(p=IncAlpha[1],mean=0,sd=1)         # compute under the null (Ho)
            thea[1] <- qnorm(p=1-IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)
        }else{
            stop("direction should be either greater or smaller")
        }
    }
    thealpha[1] <- IncAlpha[1]   
    thebeta[1] <- IncBeta[1]
    if(direction=="smaller"){
        thea <- pmin(thea,theb) # just in case of over-running
    }else{
        theb <- pmin(thea,theb) # just in case of over-running
    }
    ## if(Trace){
    ## cat("\n a.1 computed as",thea[1],"and b.1 as",theb[1]," \n")
    ## }
    ## }}}
    ## {{{ loop over k >=2
    if(K>1){
        for(k in 2:K){
            if(!thea[k-1]==theb[k-1]){
                ## {{{ if  over-running has not occurred yet
                thealpha[k] <- f(I=(k/K)*Imax,rho=rho,alpha=alpha,Imax=Imax) 
                IncAlpha[k] <- thealpha[k] - thealpha[(k-1)]   
                thebeta[k] <- g(I=(k/K)*Imax,rho=rho,beta=beta,Imax=Imax)  
                IncBeta[k] <- thebeta[k] - thebeta[(k-1)]   
                if(binding){
                    ## {{{ 
                    if(direction=="smaller"){
                        ## {{{ b_k by solving what follows 
                        theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows (when binding = TRUE )
                        try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(thea[1:(k-1)],x),
                                                                   upper = c(theb[1:(k-1)],Inf),
                                                                   mean=rep(0,k),
                                                                   sigma= sigmaZk[1:k,1:k],
                                                                   abseps = abseps) - IncAlpha[k]},
                                               lower = thea[k-1],
                                               upper = theb[k-1],
                                               tol = abseps)$root, silent = TRUE)
                        IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
                        ## }}}
                        ## {{{ a_k by solving what follows 
                        if(IsbkOK){
                            thea[k] <- theb[k] # just to handle cases in which there is no root in what follows 
                            try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(thea[1:(k-1)],-Inf),
                                                                       upper = c(theb[1:(k-1)],x),
                                                                       mean=thetheta[1:k],
                                                                       sigma= sigmaZk[1:k,1:k],
                                                                       abseps = abseps) - IncBeta[k]},
                                                   lower = thea[k-1], 
                                                   upper = theb[k], 
                                                   tol = abseps)$root, silent = TRUE)
                        }else{
                            thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
                        }
                        ## }}}
                    }
                    if(direction=="greater"){
                        ## {{{ b_k by solving what follows 
                        theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows  (when binding = TRUE )
                        try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],-Inf),
                                                                   upper = c(thea[1:(k-1)],x),
                                                                   mean=rep(0,k), 
                                                                   sigma= sigmaZk[1:k,1:k],
                                                                   abseps = abseps) - IncAlpha[k]},
                                               lower = theb[k-1],
                                               upper = thea[k-1],
                                               tol = abseps)$root, silent = TRUE)
                        IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
                        ## }}}
                        ## {{{ a_k by solving what follows 
                        if(IsbkOK){
                            thea[k] <- theb[k] # just to handle cases in which there is no root in what follows                            
                            try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],x),
                                                                       upper = c(thea[1:(k-1)],Inf),
                                                                       mean=thetheta[1:k],
                                                                       sigma= sigmaZk[1:k,1:k],
                                                                       abseps = abseps) - IncBeta[k]},
                                                   lower = theb[k],
                                                   upper = thea[k-1],
                                                   tol = abseps)$root, silent = TRUE)
                        }else{
                            thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
                        }
                        ## }}}
                    }
                    ## }}}
                }else{
                    ## {{{ 
                    if(direction=="smaller"){
                        ## {{{ b_k by solving what follows
                        theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows
                        try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(rep(-Inf, k-1),x),
                                                                   upper = c(theb[1:(k-1)],Inf),                                                                   
                                                                   mean=rep(0,k),
                                                                   sigma= sigmaZk[1:k,1:k],
                                                                   abseps = abseps) - IncAlpha[k]},
                                               lower = thea[1],
                                               upper = theb[1],
                                               tol = abseps)$root, silent = TRUE)
                        ## browser()
                        IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
                        ## }}}
                        ## {{{ a_k by solving what follows
                        if(IsbkOK){
                            thea[k] <- theb[k] # just to handle cases in which there is no root in what follows 
                            try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(thea[1:(k-1)],-Inf),
                                                                       upper = c(theb[1:(k-1)],x),
                                                                       mean=thetheta[1:k],
                                                                       sigma= sigmaZk[1:k,1:k],
                                                                       abseps = abseps) - IncBeta[k]},
                                                   lower = thea[1],
                                                   upper = theb[1],
                                                   tol = abseps)$root, silent = TRUE)
                        }else{
                            thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
                        }
                        ## }}}
                    }
                    if(direction=="greater"){
                        ## {{{ b_k by solving what follows
                        ## browser()
                        theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows
                        try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],-Inf),
                                                                   upper = c(rep(Inf, k-1),x),
                                                                   mean=rep(0,k), #thetheta[1:k],
                                                                   sigma= sigmaZk[1:k,1:k],
                                                                   abseps = abseps) - IncAlpha[k]},
                                               lower = theb[k-1],
                                               upper = thea[k-1],
                                               tol = abseps)$root, silent = TRUE)
                        IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
                        ## }}}
                        ## {{{ a_k by solving what follows
                        if(IsbkOK){
                            thea[k] <- theb[k] # just to handle cases in which there is no root in what follows                            
                            try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],x),
                                                                       upper = c(thea[1:(k-1)],Inf),
                                                                       mean=thetheta[1:k],
                                                                       sigma= sigmaZk[1:k,1:k],
                                                                       abseps = abseps) - IncBeta[k]},
                                                   lower = theb[1],
                                                   upper = thea[1],
                                                   tol = abseps)$root, silent = TRUE)
                        }else{
                            thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
                        }
                        ## }}}
                    }
                    ## }}}
                }
                ## {{{ to deal with over-runnig (see chapter Jennison) 
                if(direction=="smaller"){
                    thea <- pmin(thea,theb)
                }else{
                    theb <- pmin(thea,theb)
                }
                ## }}}
                ## }}}
            }else{
                ## {{{ if  over-running has already occurred
                thea[k:K] <- thea[k-1]
                theb[k:K] <- theb[k-1]
                ## }}}
            }
        }
    }
    ## }}}
    ## {{{ create  output
    d <- data.frame(a.k=thea,
                    b.k=theb,
                    Type.I.Error=thealpha,
                    Type.II.Error=thebeta,
                    Inc.Type.I=IncAlpha,
                    Inc.Type.II=IncBeta,
                    I.k=((1:K)/K)*Imax
                    )
    out <- list(boundaries=d,
                rho=rho,
                alpha=alpha,
                beta=beta,
                K=K,
                If=If,
                Imax=Imax,
                theta=theta,
                coef=mycoef,
                abseps=abseps,
                toldiff=toldiff,
                binding=binding,
                direction=direction)
    class(out) <- "SeqCR"
    ## }}}
    .Random.seed <<- old # restore the current seed (before the call to the function)
    out
}
