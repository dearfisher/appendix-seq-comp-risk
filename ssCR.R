
ssCR <- function(NMC,          # number of Monte Carlo simulations
                 tau1,         # length of accrual period
                 tau2,         # time of end of study (i.e. maximum follow-up time)
                 maxN=100,     # maximum sample size
                 minN=50,      # minimum sample size
                 byN=10,       # increment, such that we compute for sample sizes in seq(from=minN,to=maxN,by=byN)
                 mypower=0.8,  # power 
                 alpha=0.05,   # type-I error
                 myseed=20180616, # seed, for reproducible results
                 times,           # times at which we specify the expected values of the cumulative incidence functions (CIF).
                 pTreat=0.5,      # randomization ratio, e.g. 0.5 for 1:1.
                 CIF10,        # expected values of the CIF of main event in control group, at each time.
                 CIF20,        # same, but for the competing event
                 CIF11,        # same, but for the main event in the treatment group.
                 CIF21         # same, but for the competing event
                 ){
    StartClock <- Sys.time()
    set.seed(myseed)   

    ## {{{ stop if incorrect input values
    if(max(CIF11+CIF21,CIF10+CIF20)>1){stop("Incorrect input: the sum of the cumulative icidence functions of the main and competing event cannot be larger than 1.")}
    if(any(c(length(CIF10),length(CIF20),length(CIF11),length(CIF21))!=length(times))){
        stop("Incorrect input: the length of vectors CIF10, CIF20, CIF11, CIF21 and times should be the same.")}
    if(maxN<minN){
        stop("Incorrect input: maxN<minN.")}
    ## }}}

    ## {{{ load packages
    require(prodlim)
    require(survival)
    require(cmprsk)
    ## }}})
    
    ## {{{ helper functions
    # generate time to event fom u~U([0,1]) and risk function, specified from sum of CIFs and corresponding time points
    GenTime <- function(u,allRisk,timepoints){
        i <- sum(u>=allRisk)
        if(i==0){
            out <- u*timepoints[1]/allRisk[1]
        }else{
            if(u<=max(allRisk)){
                out <- (u-allRisk[i])*(timepoints[i+1]-timepoints[i])/(allRisk[i+1]-allRisk[i]) + timepoints[i]
            }else{
                out <- max(timepoints)+1
            }
        }  
        out
    }
    # to compute survival probabilities
    fSurv <- function(t,
                      times,
                      CIF1,
                      CIF2
                      ){
        i <- sum(t > times)
        if(i==0){
            out <- 1-(t*CIF1[1]/times[1] + t*CIF2[1]/times[1])
        }else{
            if(t<=max(times)){
                out <- 1 - ( (t-times[i])*(CIF1[i+1]-CIF1[i])/(times[i+1]-times[i]) + (t-times[i])*(CIF2[i+1]-CIF2[i])/(times[i+1]-times[i])  +  (CIF1[i] + CIF2[i]))
            }else{
                out <- NA
            }
        }  
        out
    }
    # to compute the subdistribution hazard of event 1 at any time
    dsubd1 <- function(t,
                       times,
                       CIF1
                       ){
        i <- sum(t > times)
        if(i==0){
            si <- CIF1[1]/times[1] # slope, i.e. derivative
            CIFt <- t*si           # CIF at t
            out <- si/(1-CIFt)
        }else{
            if(t<=max(times)){
                si <- (CIF1[i+1]-CIF1[i])/(times[i+1]-times[i])
                CIFt <-  (t-times[i])*si + CIF1[i]
                out <- si/(1-CIFt)
            }else{
                out <- NA
            }
        }  
        out
    }
    # to compute cause specific hazards
    fCSH <- function(t,
                     times,
                     CIF1,
                     CIF2
                     ){
        i <- sum(t > times)
        if(i==0){
            out <- (CIF1[1]/times[1])/fSurv(t=t,times=times,CIF1=CIF1,CIF2=CIF2)
        }else{
            if(t<=max(times)){
                out <- ((CIF1[i+1]-CIF1[i])/(times[i+1]-times[i]))/fSurv(t=t,times=times,CIF1=CIF1,CIF2=CIF2)
            }else{
                out <- NA
            }
        }  
        out
    }
    ## }}}

    ## {{{  Preliminaries
    AllN <- seq(from=minN,to=maxN,by=byN)
    FGpval <- matrix(NA,nrow=NMC,ncol=length(AllN))
    colnames(FGpval) <- AllN
    CSHpval <- FGpval
    ## }}}

    #-- Monte-Carlo simulations
    for(simu in 1:NMC){

        ## {{{ generate data
        print(paste0("Simu ",simu," out of ",NMC))
        #
        Treat <- c(rep(0,round(maxN/2)), rep(1,maxN-round(maxN/2)))[sample(1:maxN)] 
        entryTimes <- runif(maxN,min=0,max=tau1)
        #
        ## {{{ Generate time to event
        event <- rep(NA,maxN)
        T <- rep(NA,maxN)
        U <- runif(maxN)
        for(i in which(Treat==0)){
            T[i] <- GenTime(u=U[i],allRisk=(CIF10+CIF20),timepoints=times)
        }
        for(i in which(Treat==1)){
            T[i] <- GenTime(u=U[i],allRisk=(CIF11+CIF21),timepoints=times)
        }
        for(i in 1:maxN){
            j <- sum(T[i]>=times)+1
            if(j<length(times)){
                if(Treat[i]==1){
                    ratioi <- (CIF11[j+1]-CIF11[j])/(CIF11[j+1]-CIF11[j] + CIF21[j+1]-CIF21[j] )
                }else{
                    ratioi <- (CIF10[j+1]-CIF10[j])/(CIF10[j+1]-CIF10[j] + CIF20[j+1]-CIF20[j] )
                }
            }else{
                ratioi <- 0.5
            }
            event[i] <- (2-rbinom(n=1,size=1,prob=ratioi))
        }
        ## }}}

        ## {{{ Add Censoring
        Ttilde <- pmin(tau2-entryTimes,T)
        Delta <- as.numeric(T<=tau2-entryTimes)
        Cause <- event*Delta
        head(cbind(U,T,Treat,event,Ttilde,Delta,Cause,tau2-entryTimes))
        d <- data.frame(U,T,Treat,uncensEvent=event,Delta,C=tau2-entryTimes,time=Ttilde,event=Cause)
        ## }}}

        ## }}}

        ## {{{ prelminary computation
        allt <- seq(from=0,to=max(times),length.out=100)
        allCSH10 <- sapply(allt,fCSH,times=times,CIF1=CIF10,CIF2=CIF20)
        allCSH11 <- sapply(allt,fCSH,times=times,CIF1=CIF11,CIF2=CIF21)
        allCSH20 <- sapply(allt,fCSH,times=times,CIF1=CIF20,CIF2=CIF10)
        allCSH21 <- sapply(allt,fCSH,times=times,CIF1=CIF21,CIF2=CIF11)
        allsubd10 <- sapply(allt,dsubd1,times=times,CIF1=CIF10)
        allsubd11 <- sapply(allt,dsubd1,times=times,CIF1=CIF11)            
        ## }}}
               
        ## {{{ compute p-values
        for(myn in AllN){
            dn <- d[1:myn,]
            xCSH <- summary(coxph(Surv(time,event==1)~Treat,data=dn))$coef
            CSHpval[simu,which(myn==AllN)] <- 1-pnorm(xCSH[,"z"]) # p-value for one sided test
            xFG <- crr(dn$time,
                       dn$event,
                       cov1=model.matrix(~dn$Treat)[,-1],
                       failcode=1,cencode=0)
            FGpval[simu,which(myn==AllN)] <- (1 - pnorm(abs(xFG$coef)/sqrt(diag(xFG$var))))         
        }
        ## }}}
    }
  
    ## {{{ summarize results on power
    # CSH
    CSHpower <- apply(CSHpval,2,function(x){mean(x<=alpha)})
    uCSHpower <- apply(CSHpval,2,function(x){binom.test(x=sum(x<=alpha),n=NMC)$conf.int[2]})
    lCSHpower <- apply(CSHpval,2,function(x){binom.test(x=sum(x<=alpha),n=NMC)$conf.int[1]})
    # Gray
    FGpower <- apply(FGpval,2,function(x){mean(x<=alpha)})
    uFGpower <- apply(FGpval,2,function(x){binom.test(x=sum(x<=alpha),n=NMC)$conf.int[2]})
    lFGpower <- apply(FGpval,2,function(x){binom.test(x=sum(x<=alpha),n=NMC)$conf.int[1]})   
    #--- CSH
    if(min(uCSHpower)<mypower & max(uCSHpower)>mypower){
        nup <- (mypower-max(uCSHpower[which(uCSHpower<=mypower)]))*byN/(min(uCSHpower[which(uCSHpower>mypower)]) -max(uCSHpower[which(uCSHpower<=mypower)])) + max(AllN[which(uCSHpower<=mypower)])
    }else{nup <- NA}
    if(min(lCSHpower)<mypower& max(lCSHpower)>mypower){    
        nlow <- (mypower-max(lCSHpower[which(lCSHpower<=mypower)]))*byN/(min(lCSHpower[which(lCSHpower>mypower)]) -max(lCSHpower[which(lCSHpower<=mypower)])) + max(AllN[which(lCSHpower<=mypower)])
    }else{nlow <- NA}
    if(min(CSHpower)<mypower & max(CSHpower)>mypower){        
        nest <- (mypower-max(CSHpower[which(CSHpower<=mypower)]))*byN/(min(CSHpower[which(CSHpower>mypower)]) -max(CSHpower[which(CSHpower<=mypower)])) + max(AllN[which(CSHpower<=mypower)])
    }else{nest <- NA}
    nCSH <- c(nest,nup,nlow)
    names(nCSH) <- c("Est.","lower","upper")
    #
    #--- FG
    if(min(uFGpower)<mypower & max(uFGpower)>mypower){
        nup <- (mypower-max(uFGpower[which(uFGpower<=mypower)]))*byN/(min(uFGpower[which(uFGpower>mypower)]) -max(uFGpower[which(uFGpower<=mypower)])) + max(AllN[which(uFGpower<=mypower)])
    }else{nup <- NA}
    if(min(lFGpower)<mypower & max(lFGpower)>mypower){    
        nlow <- (mypower-max(lFGpower[which(lFGpower<=mypower)]))*byN/(min(lFGpower[which(lFGpower>mypower)]) -max(lFGpower[which(lFGpower<=mypower)])) + max(AllN[which(lFGpower<=mypower)])
    }else{nlow <- NA}
    if(min(FGpower)<mypower & max(FGpower)>mypower){
        nest <- (mypower-max(FGpower[which(FGpower<=mypower)]))*byN/(min(FGpower[which(FGpower>mypower)]) -max(FGpower[which(FGpower<=mypower)])) + max(AllN[which(FGpower<=mypower)])
    }else{nest <- NA}
    nFG <- c(nest,nup,nlow)
    names(nFG) <- c("Est.","lower","upper")   
    ## }}}
    
    StopClock <- Sys.time()

    ## {{{ create output
    out <- list(CompTime=difftime(StopClock,StartClock),
                input=list(
                    NMC=NMC, 
                    tau1=tau1,
                    tau2=tau2,
                    maxN=maxN,
                    minN=minN,
                    byN=byN,  
                    mypower=mypower, 
                    alpha=alpha, 
                    myseed=myseed,
                    times=times, 
                    pTreat=pTreat,
                    CIF10=CIF10, 
                    CIF20=CIF20, 
                    CIF11=CIF11,
                    CIF21=CIF21
                ),
                res=list(CSH=list(power.est=CSHpower,
                                  power.lower=lCSHpower,
                                  power.upper=uCSHpower,
                                  n=nCSH,
                                  pvalue=CSHpval),
                         FG=list(power.est=FGpower,
                                 power.lower=lFGpower,
                                 power.upper=uFGpower,
                                 n=nFG,
                                 pvalue=FGpval),
                         ForPlot=list(allt=allt, 
                                      allCSH10=allCSH10, 
                                      allCSH11=allCSH11, 
                                      allCSH20=allCSH20, 
                                      allCSH21=allCSH21, 
                                      allsubd10=allsubd10,
                                      allsubd11=allsubd11),
                         OneDataSet=dn
                         )
                )
    class(out) <- "ssCR"
    ## }}}

    invisible(out)
}


plot.ssCR <- function(x,
                      type="everything",
                      mycol= c("forestgreen","red"),
                      AddEst=FALSE,
                      myxlim)
{
    ## {{{ preliminaries    
    if(! type %in% c("everything","Gray","CSH","data generation")){stop("Incorrect input value for type")}   
    if(type=="Gray" | type=="CSH"){
        par(mfrow=c(1,1),mai=c(0.8,0.8,0.6,0.2))       
    }    
    if(type=="data generation"){
        par(mfrow=c(3,2),mai=c(0.6,0.6,0.15,0.6))
    }
    if(type=="everything"){
        par(mfrow=c(4,2),mai=c(0.6,0.6,0.15,0.6))
    }
    ## }}}
    ## {{{ plot data generation    
    if(type=="data generation" | type=="everything"){
        ## {{{ extract values to plot
        tau1 <- x$input$tau1
        tau2 <- x$input$tau2
        times <- x$input$times
        times2 <- times[-length(times)]
        CIF10 <- x$input$CIF10
        CIF20 <- x$input$CIF20
        CIF11 <- x$input$CIF11
        CIF21 <- x$input$CIF21
        allt <- x$res$ForPlot$allt
        allCSH10 <- x$res$ForPlot$allCSH10
        allCSH11 <- x$res$ForPlot$allCSH11
        allCSH20 <- x$res$ForPlot$allCSH20
        allCSH21 <- x$res$ForPlot$allCSH21
        allsubd10 <- x$res$ForPlot$allsubd10
        allsubd11 <- x$res$ForPlot$allsubd11
        if(missing(myxlim)){myxlim <- quantile(times,probs=c(0,0.8))}
        ## }}}
        ## {{{ Compute CIF estimates based on one simulated data
        if(AddEst){
            times2 <- times[-length(times)]
            fit <- prodlim(Hist(time,event)~Treat,data=x$res$OneDataSet)
            ObsCIF10 <- predict(fit,cause=1,times=times2,newdata=data.frame(Treat=0))
            ObsCIF11 <- predict(fit,cause=1,times=times2,newdata=data.frame(Treat=1))
            ObsCIF20 <- predict(fit,cause=2,times=times2,newdata=data.frame(Treat=0))
            ObsCIF21 <- predict(fit,cause=2,times=times2,newdata=data.frame(Treat=1))
        }
        ## }}}
        ## {{{ plot CIF 
        #-- main event        
        plot(c(0,times),c(0,CIF10),col=mycol[1],
             type="l",
             axes=FALSE,
             ylim=c(0,1),
             xlim=myxlim,
             xlab="follow-up time", ylab="Cumulative incidence function",main="Main event",
             lwd=2)
        axis(1)
        axis(2,las=2)
        lines(c(0,times),c(0,CIF11),col=mycol[2],
              type="l",
              lwd=2)
        # add estimate for one data set of size maxN
        if(AddEst){
            points(times2,ObsCIF10,cex=2,col=mycol[1])
            points(times2,ObsCIF11,cex=2,col=mycol[2])
            plot(fit,cause=1,col=mycol,add=TRUE,lty=2)
        }
        legend("topleft",col=mycol,legend=c("control","active"),bty="n",lwd=2)
        #-- competing event
        plot(c(0,times),c(0,CIF20),col=mycol[1],
             type="l",
             axes=FALSE,ylim=c(0,1),
             xlim=myxlim,
             xlab="follow-up time", ylab="Cumulative incidence function",main="Competing event",
             lwd=2)
        axis(1)
        axis(2,las=2)
        lines(c(0,times),c(0,CIF21),col=mycol[2],
              type="l",
              lwd=2)
        if(AddEst){
            points(times2,ObsCIF20,cex=2,col=mycol[1])
            points(times2,ObsCIF21,cex=2,col=mycol[2])
            plot(fit,cause=2,col=mycol,add=TRUE,lty=2)
        }
        ## }}}
        ## {{{ plot CSH
        #--cause specific hazard
        ylimCSH <- range(c(allCSH10[allt <=max(myxlim)],allCSH11[allt <=max(myxlim)],allCSH20[allt <=max(myxlim)],allCSH21[allt <=max(myxlim)]))
        ylimFG <- range(c(allsubd10[allt <=max(myxlim)],allsubd11[allt <=max(myxlim)]))
        ylimCSHHR1 <- c(floor(min((allCSH11/allCSH10)[allt <=max(myxlim)])),
                        ceiling(max((allCSH11/allCSH10)[allt <=max(myxlim)]))
                        )
        ylimCSHHR2 <- c(floor(min((allCSH21/allCSH20)[allt <=max(myxlim)])),
                        ceiling(max((allCSH21/allCSH20)[allt <=max(myxlim)]))
                        )
        ylimSHR <- c(floor(min((allsubd11/allsubd10)[allt <=max(myxlim)])),
                     ceiling(max((allsubd11/allsubd10)[allt <=max(myxlim)])))
        # main event
        plot(allt,allCSH10,type="l",col=mycol[1],lwd=2,ylim=ylimCSH,
             xlim=myxlim, main="Main event",
             xlab="follow-up time", ylab="Hazard function",axes=FALSE)
        axis(1)
        axis(2,las=2)
        lines(allt,allCSH11,type="l",col=mycol[2],lwd=2)
        par(new=TRUE)
        plot(allt,allCSH11/allCSH10,type="l",axes=FALSE,
             ylim=ylimCSHHR1,
             xlim=myxlim,
             xlab="",ylab="")
        mtext(side=4,"HR",at=mean(ylimCSHHR1),line=2.5)
        axis(4,las=2)           
        # competing  event
        plot(allt,allCSH20,type="l",col=mycol[1],lwd=2,
             xlim=myxlim, main="Competing event",
             ylim=ylimCSH,
             xlab="follow-up time", ylab="Hazard function",axes=FALSE)
        lines(allt,allCSH21,type="l",col=mycol[2],lwd=2)
        axis(1)
        axis(2,las=2)
        par(new=TRUE)
        plot(allt,allCSH21/allCSH20,type="l",axes=FALSE,
             ylim=ylimCSHHR2,
             xlim=myxlim,
             xlab="",ylab="")
        mtext(side=4,"HR",at=mean(ylimCSHHR2),line=2.5)
        axis(4,las=2)
        ## }}}
        ## {{{ plot subdistribution hazard of main event               
        plot(allt,allsubd10,type="l",col=mycol[1],lwd=2,
             ylim=ylimFG,
             xlim=myxlim, main="Main event",
             xlab="follow-up time", ylab="Subdistribution Hazard function",axes=FALSE)
        axis(1)
        axis(2,las=2)
        lines(allt,allsubd11,type="l",col=mycol[2],lwd=2)
        par(new=TRUE)
        plot(allt,allsubd11/allsubd10,type="l",axes=FALSE,
             ylim=ylimSHR,
             xlim=myxlim,
             xlab="",ylab="")
        mtext(side=4,"SHR",at=mean(ylimSHR),line=2.5)
        axis(4,las=2)
        ## }}}
        ## {{{ plot censoring distribution
        plot(c(0,tau2-tau1,tau2),c(1,1,0),type="l",
             lwd=2,
             xlim=myxlim,
             main="Censoring",
             xlab="follow-up time t",
             ylab="P(C>t)",axes=FALSE)
        axis(1) 
        axis(2,las=2)
        if(tau2==Inf){abline(h=1,col="blue",lty=2,lwd=2)}

        if(AddEst){
            if(all(x$res$OneDataSet$event!=0)){
                abline(h=1,col="blue",lty=1,lwd=1)
            }else{
                fitC <- prodlim(Hist(time,event!=0)~1,reverse=TRUE,data=x$res$OneDataSet)
                plot(fitC,add=TRUE,col="blue",lty=2,lwd=2,type="surv")
            }
            abline(v=tau2-tau1,lty=2)
        }
        ## }}}
    }
    ## }}}
    ## {{{ plot CSH results
    AllN <- seq(from=x$input$minN,to=x$input$maxN,by=x$input$byN)   
    if(type=="CSH" | type=="everything"){
        # plot CSH sample size results
        plot(AllN,x$res$CSH$power.est*100,
             type="b",pch=19,ylim=c(0,100),ylab="power (%)",xlab="sample size",axes=FALSE,main="Logrank test")
        axis(1)
        axis(1,at=range(AllN))
        axis(2,las=2)
        abline(h=x$input$mypower*100,lwd=2,lty=2,col="red")
        polygon(x=c(AllN,rev(AllN)),
                y=c(x$res$CSH$power.upper,rev(x$res$CSH$power.lower))*100,col=rgb(0.5,0.5,0.5,0.5),border=NA)
        # display sample size
        if(!is.na(x$res$CSH$n["upper"])){
            segments(x0=x$res$CSH$n["upper"],x1=x$res$CSH$n["upper"],y0=-10,y1=x$input$mypower*100,lty=2)
            text(x=x$res$CSH$n["upper"],y=5,labels=floor(x$res$CSH$n["upper"]),pos=2)
        }
        if(!is.na(x$res$CSH$n["lower"])){
            segments(x0=x$res$CSH$n["lower"],x1=x$res$CSH$n["lower"],
                     y0=-10,y1=x$input$mypower*100,lty=2)
            text(x=x$res$CSH$n["lower"],y=5,labels=ceiling(x$res$CSH$n["lower"]),pos=4)
        }
        if(!is.na(x$res$CSH$n["Est."])){
            segments(x0=x$res$CSH$n["Est."],x1=x$res$CSH$n["Est."],y0=-10,y1=x$input$mypower*100,lty=2,col="red")
            text(x=x$res$CSH$n["Est."],y=10,labels=ceiling(x$res$CSH$n["Est."]),pos=4,col="red")
        }        
    }
    ## }}}
    ## {{{ plot FG results
    if(type=="Gray" | type=="everything"){
        # plot FG sample size results
        plot(AllN,x$res$FG$power.est*100,
             type="b",pch=19,ylim=c(0,100),ylab="power (%)",xlab="sample size",axes=FALSE,main="Gray's test")
        axis(1)
        axis(1,at=range(AllN))
        axis(2,las=2)
        abline(h=x$input$mypower*100,lwd=2,lty=2,col="red")
        polygon(x=c(AllN,rev(AllN)),
                y=c(x$res$FG$power.upper,rev(x$res$FG$power.lower))*100,col=rgb(0.5,0.5,0.5,0.5),border=NA)
        # display sample size
        if(!is.na(x$res$FG$n["upper"])){
            segments(x0=x$res$FG$n["upper"],x1=x$res$FG$n["upper"],y0=-10,y1=x$input$mypower*100,lty=2)
            text(x=x$res$FG$n["upper"],y=5,labels=floor(x$res$FG$n["upper"]),pos=2)
        }
        if(!is.na(x$res$FG$n["lower"])){
            segments(x0=x$res$FG$n["lower"],x1=x$res$FG$n["lower"],
                     y0=-10,y1=x$input$mypower*100,lty=2)
            text(x=x$res$FG$n["lower"],y=5,labels=ceiling(x$res$FG$n["lower"]),pos=4)
        }
        if(!is.na(x$res$FG$n["Est."])){
            segments(x0=x$res$FG$n["Est."],x1=x$res$FG$n["Est."],y0=-10,y1=x$input$mypower*100,lty=2,col="red")
            text(x=x$res$FG$n["Est."],y=10,labels=ceiling(x$res$FG$n["Est."]),pos=4,col="red")
        }        
    }
    ## }}}   
}

