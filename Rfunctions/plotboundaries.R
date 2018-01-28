plotboundaries <- function(SeqCR.object,
                           myylim,
                           add=FALSE,
                           col="red",
                           colObs="red",
                           xaxis=TRUE,
                           yaxis=TRUE,
                           ylagText=c(0,0,0),  # to manually adjust the location of (i.e. move up or down) the texts "Continue", "Reject" and "Accept"
                           mypch=19){ 
    ## browser()
    d <- SeqCR.object$boundaries
    ## {{{ Plot planned boundaries    
    # If planned boundaries only
    if(!("Obs.I.k" %in% names(d)) ){
        plotPlanned(SeqCR.object=SeqCR.object,
                    myylim=myylim,
                    add=add,
                    col=col,
                    lty=1,
                    yaxis=yaxis,
                    xaxis=xaxis,
                    ylagText=ylagText,
                    mypch=mypch)
    }
    # If we have some results to plot
    if("Obs.I.k" %in% names(d) ){
        if(missing(myylim)){
            myylim <- 1.1*c(min(c(d$a.k,d$Obs.a.k,d$Obs.z),na.rm=TRUE),max(c(d$b.k,d$Obs.b.k,d$Obs.z),na.rm=TRUE))
        }
        ## browser()
        plotPlanned(SeqCR.object=SeqCR.object,
                    myylim=myylim,
                    add=add,
                    col=col,
                    transparency=c(0.05,0.1),
                    lwd=1,
                    addlegend=FALSE,
                    lty=2,
                    yaxis=yaxis,
                    xaxis=xaxis,
                    ylagText=ylagText,
                    mypch=mypch)
        ## browser()
        ## {{{ Add results to the plot
        plotAnalyses(SeqCR.object,
                     myylim=myylim,
                     transparency=c(0.25,0.125),
                     col=colObs,
                     ylagText=ylagText)
        ## }}}
    }
    ## }}}
}
