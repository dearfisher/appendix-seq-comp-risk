print.SeqCR <- function(SeqCR.object, digits=3){
    ## browser()
    d <- SeqCR.object$boundaries
    alpha <- SeqCR.object$alpha
    beta <- SeqCR.object$beta
    K <- SeqCR.object$K
    rho <- SeqCR.object$rho
    theta <- round(SeqCR.object$theta,digits)
    If <- round(SeqCR.object$If,digits)
    Imax <- round(SeqCR.object$Imax,digits)
    coef <- round(SeqCR.object$coef,digits)
    Binding <- ifelse(SeqCR.object$binding,
                      "Binding futility rules.\n", "Non Binding futility rules.\n")
    Direction <- ifelse(SeqCR.object$direction=="smaller",
                        "\nOne-sided testing approach with H0: theta < 0.\n",
                        "\nOne-sided testing approach with H0: theta > 0.\n")
        cat(paste("\nBoundaries",
              " for planning a group sequential clinical trial with K=",K,
              " analyses. \nType-I error=",alpha,
              ", power=", 1-beta," (i.e. Type-II error=",beta ,
              "), rho=",rho,", theta=",theta,",\n",
              "If=",If,", Imax=",Imax,", coef=Imax/If=",coef,".\n",
              sep=""))
    cat("\n")
    cat(paste(Binding, Direction, sep=""))
    cat("\n")
    cat("(a.k= futility boundaries, b.k= efficacy boundaries)")    
    
    if("Obs.I.k" %in% names(d)){
        ## {{{ If we have already some results
        cat("\n \nInterim results and updated boundaries are shown. \n")
        # Results of the k-th analysis are included
        ## }}}
    }
    cat("\n")
    cat("\n")
    if(Direction=="smaller"){
        print(d[,1:ncol(d)], digits=digits+1)
    }else{
        print(d[,c(2:1,3:ncol(d))], digits=digits+1)
    }
}
