plotPlanned <- function(SeqCR.object,
                        myylim,
                        add=FALSE,
                        col="red",
                        transparency=c(0.25,0.125),
                        lwd=2,
                        addlegend=TRUE,
                        addvline=TRUE,
                        xaxis=TRUE,
                        yaxis=TRUE,
                        lty=1,
                        ylagText=c(0,0,0),
                        mypch=19){
    ## browser()
    ## {{{ get information to plot
    d <- SeqCR.object$boundaries
    K <- nrow(d)
    direc <- SeqCR.object$direction
    ## }}}
    ## {{{ Find appropriate ylim
    if(missing(myylim)){
        myylim <- 1.1*c(min(d$a.k),max(d$b.k))
    }
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
    ## {{{ Make plot
    if(!add){
        ## {{{ to have appropriate xlim
        if(!("Obs.I.k" %in% names(d)) ){
            xlimmax <- ceiling(max(d$I.k))
        }else{
            ## browser()
            xlimmax <- ceiling(max(c(d$I.k,d$Obs.I.k), na.rm=TRUE))
        }
        ## }}}
        ## browser()
        plot(d$I.k,d$a.k,
             type="b",xlab="Information",
             ylab="z",
             ylim=myylim,
             xlim=c(0,xlimmax),
             axes=FALSE,col=col,lwd=lwd,lty=lty,
             pch=mypch)
        if(xaxis){axis(1,at=round(seq(from=0,to=ceiling(max(d$I.k)),length.out=K+1),2))}
        if(yaxis){axis(2,las=2)}
    }else{
        lines(d$I.k,d$a.k,type="b",col=col,lwd=lwd,lty=lty,pch=mypch)
    }
    lines(d$I.k,d$b.k,type="b",col=col,lwd=lwd,lty=lty,pch=mypch)
    abline(h=0,lty=2,lwd=1)
    if(addvline){
        segments(y0=d$a.k,y1=rep(alim,K),x0=d$I.k,x1=d$I.k,lty=3,col=col,lwd=lwd)
        segments(y0=d$b.k,y1=rep(blim,K),x0=d$I.k,x1=d$I.k,lty=3,col=col,lwd=lwd)
    }
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
    polygon(y=c(d$a.k,rep(alim,K)),x=c(d$I.k,rev(d$I.k)),col=myrgbcol1,border=NA)
    polygon(y=c(d$b.k,rep(blim,K)),x=c(d$I.k,rev(d$I.k)),col=myrgbcol2,border=NA)
    ## }}}
    ## {{{ Add text to plot
    if(addlegend){
        text(x=0,y=min(myylim),"Analysis:",pos=4)
        text(x=d$I.k,y=rep(min(myylim),K),pos=4,labels = 1:K)
        text(0,d$b.k[K]+ylagText[1],"Continue",pos=4)    
        text(d$I.k[K-1]+0.5*(d$I.k[K]-d$I.k[K-1]),
             d$b.k[K]+0.5*(max(myylim)-d$b.k[K])+ylagText[2],
             ifelse(SeqCR.object$direction=="smaller",
                    "Reject H0","Accept H0"))
        text(d$I.k[K-1]+0.5*(d$I.k[K]-d$I.k[K-1]),
             d$b.k[K]+0.5*(min(myylim)-d$b.k[K])+ylagText[3],
             ifelse(SeqCR.object$direction=="smaller","Accept H0","Reject H0"))
    }
    ## }}}
}
