plot.SeqCR <- function(SeqCR.object,
                       myylim=c(-4,4),
                       col="blue",
                       colObs="red"){
    # call plotboundaries fucntion
    plotboundaries(SeqCR.object,
                   myylim=myylim,
                   col=col,
                   colObs=colObs)
}
