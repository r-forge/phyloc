"analyse.DI" <-
function(phylogeny, data, log.interval=FALSE, nbins=NULL){

#Performs difference vs interval analyses for a phylogeny.

#10/5/06: Produced de novo; links to matched.pairs.DI and binned.DI

phylogeny.mat<-cophenetic(phylogeny)
data.mat<-as.matrix(dist(data))
all.matched.pairs<-find.independent.pairs(phylogeny)

di<-matched.pairs.DI(phylogeny, phylogeny.mat, all.matched.pairs, data.mat)
bdi<-binned.DI(di, log.interval=log.interval, nbins=nbins)

weighted.model<-lm(bdi[,1]~bdi[,2],weights=bdi[,3])
unweighted.model<-lm(bdi[,1]~bdi[,2])
subhead<-paste("Wt slope =",format(as.numeric(coef(weighted.model)[2])),": Unwt slope = ",format(as.numeric(coef(unweighted.model)[2])))
plot(bdi[,1]~bdi[,2],xlab=names(bdi)[2],ylab=names(bdi)[1],cex.lab=1.3,sub=subhead, main=paste("Sample size =",sum(bdi$Sample.Size)))
abline(weighted.model)
abline(unweighted.model,lty=2)

to.return<-list(weighted.model,unweighted.model)

return(to.return)

}

