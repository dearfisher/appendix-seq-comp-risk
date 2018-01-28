plotAnalyses <- function(SeqCR.object,
                         col="red",
                         transparency=c(0.1,0.2),
                         lwd=2,
                         myylim,
                         ylagText=c(0,0,0),
                         mypch=19
                         ){
    ## browser()
    ## {{{ get parameters
    d <- SeqCR.object$boundaries
    K <- nrow(d)
    direc <- SeqCR.object$direction
    ## }}}
    ## {{{ velues depending on direction
    if(direc=="smaller"){
        alim <- -20
        blim <- 20
    }else{
        alim <- +20
        blim <- -20
    }
    ## }}}
    ## {{{ Add results to the plot
    ## browser()
    #interim analysis results
    ## points(d$Obs.I.k,d$Obs.a.k,col=col,lwd=lwd,pch=mypch)
    ## points(d$Obs.I.k,d$Obs.b.k,col=col,lwd=lwd,pch=mypch)
    lines(d$Obs.I.k,d$Obs.a.k,col=col,lwd=lwd,pch=mypch,type="b")
    lines(d$Obs.I.k,d$Obs.b.k,col=col,lwd=lwd,pch=mypch,type="b")
    ## browser()
    points(d$Obs.I.k,d$Obs.z.k,col=col,lwd=lwd,
           pch=4,
           cex=2)

    colvect <- col2rgb(col)/255
    myrgbcol1 <- rgb(red=colvect[1],
                     green=colvect[2],
                     blue=colvect[3],
                     alpha=transparency[1])
    myrgbcol2 <- rgb(red=colvect[1],
                     green=colvect[2],
                     blue=colvect[3],
                     alpha=transparency[2])
    ## browser()    
    polygon(y=c(d$Obs.a.k[!is.na(d$Obs.a.k)],rep(alim,sum(!is.na(d$Obs.a.k)))),
            x=c(d$Obs.I.k[!is.na(d$Obs.a.k)],rev(d$Obs.I.k[!is.na(d$Obs.a.k)])),
            col=myrgbcol1,
            border=myrgbcol1)
    polygon(y=c(d$Obs.b.k[!is.na(d$Obs.b.k)],rep(blim,sum(!is.na(d$Obs.b.k)))),
            x=c(d$Obs.I.k[!is.na(d$Obs.b.k)],rev(d$Obs.I.k[!is.na(d$Obs.b.k)])),
            col=myrgbcol2,
            border=myrgbcol1)
    segments(y0=d$Obs.a.k[!is.na(d$Obs.a.k)],
             y1=rep(alim,sum(!is.na(d$Obs.a.k))),
             x0=d$Obs.I.k[!is.na(d$Obs.a.k)],
             x1=d$Obs.I.k[!is.na(d$Obs.a.k)],
             lty=3,col=col,lwd=lwd)
    segments(y0=d$Obs.b.k[!is.na(d$Obs.b.k)],
             y1=rep(blim,sum(!is.na(d$Obs.b.k))),
             x0=d$Obs.I.k[!is.na(d$Obs.b.k)],
             x1=d$Obs.I.k[!is.na(d$Obs.b.k)],
             lty=3,col=col,lwd=lwd)
    ## }}}
    ## {{{ Add legend
    ## browser()
    text(x=0,y=min(myylim),"Analysis:",pos=4)
    text(x=d$Obs.I.k[!is.na(d$Obs.I.k)],
         y=rep(min(myylim),sum(!is.na(d$Obs.I.k))),
         pos=4,labels = 1:sum(!is.na(d$Obs.I.k)))
    ## browser()
    text(d$Obs.I.k[1]/2,
    (d$Obs.b.k[1]+d$Obs.a.k[1])/2 + ylagText[1],
    "Continue",pos=4)
    # case where more than one interim analysis
    if(sum(!is.na(d$Obs.b.k)) > 1){
        xtext <- d$Obs.I.k[sum(!is.na(d$Obs.b.k))] - 0.5*(d$Obs.I.k[sum(!is.na(d$Obs.b.k))]-d$Obs.I.k[sum(!is.na(d$Obs.b.k))-1])
        ytext1 <- d$Obs.b.k[1]
        ytext2 <- d$Obs.a.k[1]
        text(xtext,
             ytext1 + ylagText[2],
             ifelse(direc=="smaller","Reject H0","Accept H0"))    
        text(xtext,
             ytext2 + ylagText[3],
             ifelse(direc=="smaller","Accept H0","Reject H0"))
    }else{
        # case where only one interim analysis
        ## xtext <- 1.5*d$Obs.I.k[sum(!is.na(d$Obs.b.k))]
        text(d$I.k[K-1]+0.5*(d$I.k[K]-d$I.k[K-1]),
             d$b.k[K]+0.5*(max(myylim)-d$b.k[K]) + ylagText[2],
             ifelse(direc=="smaller","Reject H0","Accept H0"))
        text(d$I.k[K-1]+0.5*(d$I.k[K]-d$I.k[K-1]),
             d$b.k[K]+0.5*(min(myylim)-d$b.k[K]) + ylagText[3],
             ifelse(direc=="smaller","Accept H0","Reject H0"))
        
    }
    ## browser()
    ## }}}
}
