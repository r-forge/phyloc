"binned.DI" <-
function(di, log.interval=FALSE, nbins=NULL){

#Takes output from matched.pairs.DI, groups the observations into bins, fits a linear
#regression to the within-bin means, plots the within-bin means and line of best fit, and
#returns the regression model.

#10/5/06: Created de novo.

di<-subset(di,select=c(Difference, Interval))

if (log.interval==TRUE){
	di$Interval<-log10(di$Interval)
	names(di)[2]<-"Log10.Interval"
	}

if (is.null(nbins)) nbins<-log2(length(di[,1]))+1 #Equivalent to Sturges' formula for numbers of histogram classes
bin<-as.factor(round(0.49999+(c(1:length(di[,1]))/(length(di[,1])/nbins))))
mean.diff<-as.numeric(tapply(di[,1],bin,mean))
var.diff<-as.numeric(tapply(di[,1],bin,var)) #uses n-1 as denominator
n.diff<-as.numeric(tapply(di[,1],bin,length))
se.diff<-sqrt(var.diff)/n.diff
mean.interval<-as.numeric(tapply(di[,2],bin,mean))
wt<-1/(se.diff)^2

to.return<-data.frame(mean.diff,mean.interval,wt,n.diff)
names(to.return)<-c(names(di),"Weight","Sample.Size")

return(to.return)
}

