#Function to plot time to event data

#Input :
#
# Time: the time to event variable
# Event: event indicator (0: censoring, 1: main event, 2: competing event)
# Date: start of time to event
# Group: (treatment) grouping variable (data will be displayed with different colors per group)
# WhatToPlot: calendar time or time in study
# dateInterim: date of an interim analysis
# OnlyAvailable: plot only the time to event data that is available at the interim analysis
# cex: size (cex) of the points on the plot
# color : vector of color to display the different groups (should be in the same order as the group appear in the data)
#
#


PlotSeqTrialTTE <- function(time,
                            event,
                            date,
                            group=NA,
                            WhatToPlot,
                            dateInterim=NA, 
                            OnlyAvailable=FALSE, 
                            legend=T,
                            cex=2,
                            color=NA,
                            showdays=TRUE){
    ## browser()
    DateStart <- min(date,na.rm=T)
    DateStop <- max(date + time,na.rm=T) + 1
    ## {{{ if we do not want to show the eact date, but just year/month
    if(showdays){
        PrintDateStart <- format(as.Date(DateStart), "%Y-%m-%d")
        PrintDateStop <- format(as.Date(DateStop), "%Y-%m-%d")
    }else{
        PrintDateStart <- format(as.Date(DateStart), "%Y-%m")
        PrintDateStop <- format(as.Date(DateStop), "%Y-%m")
    }
    ## }}}
    TimeLength <- as.numeric(difftime(DateStop,DateStart,units="days"))
    MaxTime <- max(time,na.rm=T)
    # create vector of pch symbls to plot for each subject
    mypch <- rep(1,length(event)) # pch=1 for censored
    mypch[which(event==1)] <- 19 # pch=19 for main event
    mypch[which(event==2)] <- 17 # pch=17 for competing event       
    ## npatient <- length(time)
    if(!is.na(dateInterim)){
        timeDataInterim <- as.numeric(difftime(dateInterim,DateStart,units="days"))
    }
    ## browser()
    # color of lines to display
    if(missing(color) | length(color)!=length(unique(group))){
        color <- 1:length(unique(group))
    }
    colorVector <- rep(NA,length(time))
    ## browser()
    for(i in 1:length(color)){
        colorVector[which(group==unique(group)[i])] <- color[i]
    }
    ## browser()
    if(OnlyAvailable){
        if(!is.na(dateInterim)){
            ## time <- pmin(time + as.numeric(difftime(date,DateStart,units="days")),timeDataInterim) - as.numeric(difftime(date,DateStart,units="days"))
            ## browser()
            event2 <- as.numeric(time<=as.numeric(difftime(dateInterim,date,units="days")))*event
            #table(event2)
            #table(event)
            difftimeToSubstr <- pmin(time,as.numeric(difftime(dateInterim,date,units="days")))
            time2 <- pmax(difftimeToSubstr,0)
            cbind(time,time2)
            time <- time2[difftimeToSubstr>0]
            event <- event2[difftimeToSubstr>0]
            date <- date[difftimeToSubstr>0]
        }else(stop("Please provide me with something for the argument dateInterim"))
    }
    # plot all data : calendar time
    if(WhatToPlot=="AllCalendar"){
        plot(1, type="n",
             xlab="Calendar Time",
             ylab="Patient id",
             xlim=c(0, TimeLength),
             ylim=c(0,length(time)),
             axes=FALSE)
        if(is.na(dateInterim)){
            axis(1,at=c(0,TimeLength),labels=c(PrintDateStart,PrintDateStop))
        }else{
            axis(1,at=c(0,timeDataInterim,TimeLength),labels=c(PrintDateStart,dateInterim,PrintDateStart))
            abline(v=timeDataInterim,lty=2,col="red")
        }
        axis(2,at=c(1,length(time)),las=2)
        ## browser()
        for(i in 1:length(time)){
            timelag <- as.numeric(difftime(date[i],DateStart,units="days"))
            segments(x0=timelag,x1=timelag+time[i],y0=i,y=i,col=colorVector[i])
            points(x=timelag+time[i],y=i,pch=mypch[i],cex=cex,col=colorVector[i])
        }
    }
    # plot all data : time into study
    if(WhatToPlot=="AllTimeToEvent"){
        plot(1, type="n",
             xlab="Time to event (days)",
             ylab="Patient id",
             xlim=c(0,MaxTime),
             ylim=c(0,length(time)),
             axes=FALSE)
        axis(1)
        ## axis(1,at=c(0,MaxTime))
        axis(2,at=c(1,length(time)),las=2)
        ## browser()
        for(i in 1:length(time)){
            timelag <- 0
            segments(x0=timelag,x1=timelag+time[i],y0=i,y=i,col=colorVector[i])
            points(x=timelag+time[i],y=i,pch=mypch[i],cex=cex,col=colorVector[i])
        }
    }
    if(legend){
        whichInLEgend <- which(c(0,1,2) %in% event)
        legend("bottomright",
               pch=c(1,19,17)[whichInLEgend],
               c("Censored","Main event","Competing event")[whichInLEgend],
               bty="n"
               )
        if(!missing(group)){
            legend("topleft",
                   fill = color,
                   legend=unique(group),bty="n",
                   title="Group:"
                   )
        }
    }
    #return(cbind(time,event))
}
