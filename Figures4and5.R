rm(list=ls())

library(mvtnorm)

setwd("~/path/to/your/directory")

## {{{ source functions
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir("Rcode/Rfunctions/")
## }}}


## {{{ Parameters
# group sequential trial design parameters
myrho <- 1:3
myalpha <- 0.05
mybeta <- 0.2
myK <- 2
mytheta <- log(2)
mybinding <- FALSE
# graphical paramters
mydirection <- "smaller"
allcol <- c("black","grey40","grey60")
allpch <- 15:17
myplotcoef <- 0.7         # to control the size of the plot
WheretoSavePlot <- "Fig/" #  where to save the plots
SavePlot <- TRUE          # whether to save the plots
## }}}


## {{{ Compute boundaries with different values for rho
runif(1) # jut to initiate random seed
# Compute boundaries with rho=myrho[1]
MyResults1 <- planboundaries(rho=myrho[1],  
                             alpha=myalpha,
                             beta=mybeta,
                             K=myK,         
                             theta=mytheta,
                             abseps = 1e-06,
                             binding=mybinding,
                             direction=mydirection,
                             Trace=TRUE,
                             nWhileMax=30,
                             toldiff=1e-05,
                             mycoefMax= 1.5,       
                             mycoefL=1             
                             )
# Compute boundaries with rho=myrho[2]
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
# Compute boundaries with rho=myrho[3]
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
## }}}


## {{{ Plot
if(SavePlot){
    pdf(paste0(WheretoSavePlot,"illustrate-rho-choice.pdf"),
        width=myplotcoef*10,height=myplotcoef*8)
}
par(mai=c(0.85,0.85,0.85,0.2))
# Plot planned boundaries
plotPlanned(MyResults1,
            myylim=c(-0.5,3),
            transparency=c(0,0),
            addlegend=FALSE,
            addvline=FALSE,
            col=allcol[1],
            xaxis=FALSE,
            yaxis=FALSE,
            mypch=allpch[1])
plotPlanned(MyResults2,
            add=TRUE,
            col=allcol[2],
            addlegend=FALSE,
            transparency=c(0,0),
            addvline=FALSE,
            mypch=allpch[2])
plotPlanned(MyResults3,
            add=TRUE,
            col=allcol[3],
            addlegend=FALSE,
            transparency=c(0,0),
            addvline=FALSE,
            mypch=allpch[3])
# Add non-sequential design
points(x=MyResults3$If,y=qnorm(1-myalpha),pch=4,cex=1.5,lwd=2)
segments(x0=10.5,x1=50,y0=0,y1=0,col="white",lwd=2) # minor cosmetic
# Compute increase in maximal sample size (to print in legend)
IncNSubj <- round(100*ceiling(4*c(MyResults1$boundaries$I.k[myK],
                                  MyResults2$boundaries$I.k[myK],
                                  MyResults3$boundaries$I.k[myK]))/ceiling(MyResults3$If*4) - 100,1)
# Create legend text
mylengend <- c(bquote(paste(rho==1,",  ",.( format(IncNSubj[1],digits=0,nsmall=0)),"%")),
               bquote(paste(rho==2,",  ",.( format(IncNSubj[2],digits=0,nsmall=0)),"%")),
               bquote(paste(rho==3,",  ",.( format(IncNSubj[3],digits=0,nsmall=0)),"%")),
               "Non sequential")
legend("bottomright",col=c(allcol,"black"),lwd=2,lty=c(rep(1,3),NA),pch=c(allpch,4),
       legend=as.expression(mylengend),title="Increase in maximal \n number of events required",bty="n")
# Compute required number of events
myI4 <- ceiling(c(0,MyResults1$boundaries$I.k[c(1)],MyResults3$If,MyResults1$boundaries$I.k[myK])*4)
# Add axes
axis(side=3,at=myI4/4,myI4)
axis(side=2,at=c(-0.5,0,1,2,3),labels=c("-0.5","0","1","2","3"),las=2)
axis(side=1,at=round(c(0,MyResults1$boundaries$I.k[c(1)],
                       MyResults3$If,MyResults1$boundaries$I.k[myK]),1),
     round(c(0,MyResults1$boundaries$I.k[c(1)],MyResults3$If,MyResults1$boundaries$I.k[myK]),1))
mtext("Number of events (i.e. = 4 x Information)", side=3, line=2.2)

if(SavePlot){
    dev.off()
}
## }}}


## {{{ Compute average number of required events before stoping the trial

# Values of the fraction of the expected effect size for which we
# want to compute the expected number of required events before stoping the trial
myvalues <- seq(from=0,to=1.5,by=0.05)
# With rho=myrho[1]
ResAvNCont1 <- ExpectSampleSize(MyResults1,
                                actual.binding = TRUE,
                                thetaValues=myvalues)
# Remark: we assume that we never (respect. always) decide to continue even if we
# can stop for futility, when we set "actual.binding" to "FALSE" (repect. TRUE)
# With rho=myrho[2]
ResAvNCont2 <- ExpectSampleSize(MyResults2,
                                actual.binding = TRUE,
                                thetaValues=myvalues)
# With rho=myrho[3]
ResAvNCont3 <- ExpectSampleSize(MyResults3,
                                actual.binding = TRUE,
                                thetaValues=myvalues) 
# Results when assuming that the treatmnt effect is indeed theta=delta (i.e. H_1)
AllAvNEventsH1 <- c(round(ResAvNCont1$N*(MyResults1$If*4)/100,1)[which(myvalues==1)],
                 round(ResAvNCont2$N*(MyResults2$If*4)/100,1)[which(myvalues==1)],
                 round(ResAvNCont3$N*(MyResults3$If*4)/100,1)[which(myvalues==1)])
names(AllAvNEventsH1) <- paste("rhol=",1:3,sep="")
# Results when assuming that the treatmnt effect is actually theta=0 (i.e. H_0)
AllAvNEventsH0 <- c(round(ResAvNCont1$N*(MyResults1$If*4)/100,1)[which(myvalues==0)],
                    round(ResAvNCont2$N*(MyResults2$If*4)/100,1)[which(myvalues==0)],
                    round(ResAvNCont3$N*(MyResults3$If*4)/100,1)[which(myvalues==0)])
names(AllAvNEventsH0) <- paste("rhol=",1:3,sep="")

# Plot
if(SavePlot){
    pdf(paste0(WheretoSavePlot,"illustrate-Average-events.pdf"),width=myplotcoef*10,height=myplotcoef*7)
}
par(mai=c(0.85,0.85,0,0.1))
myylim <- c(32,53)
plot(myvalues,ResAvNCont1$N*(MyResults1$If*4)/100,
     col=allcol[1],type="l",lwd=2,ylim=myylim,
     ylab="Expected number of events",
     xlab=expression(paste("Fraction of the value of the expected effect size ",delta)),
     axes=FALSE
     )
# Emphasize results under H0 and H1 (when the treatment effect corresponds to the expected one)
segments(x0=0,x1=0,y0=0,y1=max(AllAvNEventsH0),col="grey60",lty=2,lwd=1)
segments(x0=1,x1=1,y0=0,y1=max(AllAvNEventsH1),col="grey60",lty=2,lwd=1)
segments(x0=-10,x1=0,y0=AllAvNEventsH0[1],y1=AllAvNEventsH0[1],col="grey60",lty=2,lwd=1)
segments(x0=-10,x1=0,y0=AllAvNEventsH0[2],y1=AllAvNEventsH0[2],col="grey60",lty=2,lwd=1)
segments(x0=-10,x1=0,y0=AllAvNEventsH0[3],y1=AllAvNEventsH0[3],col="grey60",lty=2,lwd=1)
segments(x0=-10,x1=1,y0=AllAvNEventsH1[1],y1=AllAvNEventsH1[1],col="grey60",lty=2,lwd=1)
segments(x0=-10,x1=1,y0=AllAvNEventsH1[2],y1=AllAvNEventsH1[2],col="grey60",lty=2,lwd=1)
segments(x0=-10,x1=1,y0=AllAvNEventsH1[3],y1=AllAvNEventsH1[3],col="grey60",lty=2,lwd=1)
points(x=rep(0,3),y=AllAvNEventsH0,col=allcol,lty=2,lwd=1,pch=allpch)
points(x=rep(1,3),y=AllAvNEventsH1,col=allcol,lty=2,lwd=1,pch=allpch)
lines(myvalues,ResAvNCont1$N*(MyResults1$If*4)/100,
      col=allcol[1],type="l",lwd=2)
lines(myvalues,ResAvNCont2$N*(MyResults1$If*4)/100,
      col=allcol[2],type="l",lwd=2)
lines(myvalues,ResAvNCont3$N*(MyResults1$If*4)/100,
      col=allcol[3],type="l",lwd=2)
# Add the case of non sequential trial 
abline(h=ceiling(MyResults3$If*4),lty=2,lwd=2)
# Add axes with values of key results
axis(1,at=c(0,0.5,1,1.5),labels=c("0","1/2","1","3/2"))
axis(2,at=sort(c(myylim[1],ceiling(MyResults1$If*4),
                 AllAvNEventsH0,
                 AllAvNEventsH1)),
     las=2,labels=c(format(round(myylim[1]), nsmall=0),
                    format(sort(c( AllAvNEventsH0,
                                  AllAvNEventsH1)), nsmall=1),
                    format(ceiling(MyResults1$If*4),nsmall=0)))
# Emphasize situation corresponding to H0 and H1 (when the treatment effect corresponds to the expected one)
text(0.12,myylim[1],expression(H[0] : theta==0),cex=1.2)
text(1.12,myylim[1],expression(H[1] : theta==delta),cex=1.2)
# Add legend
mylengend <- c(bquote(paste(rho==1)),
               bquote(paste(rho==2)),
               bquote(paste(rho==3)),
               "Non sequential")
legend(x=1.1,y=52,
       col=c(allcol,"black"),
       lwd=2,
       lty=c(rep(1,3),2),
       pch=c(allpch,NA),
       legend=as.expression(mylengend),bty="n"
       )
if(SavePlot){
    dev.off()
}
## }}}
