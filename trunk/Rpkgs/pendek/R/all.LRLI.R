"all.LRLI" <-
function(phylogeny, data,plot="n"){

#Produces Gingerich-type log-rate log-interval plot for a phylogeny
#Returns three lm objects as a list
#  1. Highly nonindependent analysis using all-by-all comparisons (twice!)
#  2. Independent analysis of cherries only
#  3. Independent analysis of maximum number of pairs

#10/5/06: Bug fix to prevent crash when two species have exactly the same trait value and
#so have zero rate.  See also matched.pairs.LRLI

#10/5/06: Change to reflect change in matched.pairs.LRLI, which now returns data rather than a model.

#11/5/06: Realisation that the function assumes the species data are in the same order as tip.label.
#I might want to consider either using row.names or passing a data frame instead of a vector.

phylogeny.mat<-cophenetic(phylogeny)
data.mat<-as.matrix(dist(data))
rate.mat<-data.mat/phylogeny.mat

#Non-phylogenetic analysis
#Note that, in addition to the obvious problem of nonindependence, both halves of the matrix are used
#so there are twice as many points as there should be
rate<-as.vector(rate.mat)
rate[rate==0]<-min(rate[rate>0])/10 #Hack
interval<-as.vector(phylogeny.mat)

log10.rate<-log10(rate)
log10.interval<-log10(interval)

model1<-lm(log10.rate~log10.interval)

if (plot == "y")
{
	par(mfrow=c(3,1))
	plot(log10.rate~log10.interval,ylab="Log10(Rate)",xlab="Log10(Interval)",cex.lab=1.3)
	abline(model1)
}

#Phylogenetic analysis of cherries only
cherries<-find.cherries(phylogeny)
ds<-matched.pairs.LRLI(phylogeny, phylogeny.mat, cherries, data.mat)
model2<-lm(ds$Log10.Rate~ds$Log10.Interval)

if (plot == "y")
{
	plot(ds$Log10.Rate~ds$Log10.Interval,ylab="Log10(Rate)",xlab="Log10(Interval)",cex.lab=1.3)
	abline(model2)
}

#Phylogenetic analysis of all independent pairs
all.matched.pairs<-find.independent.pairs(phylogeny)
ds<-matched.pairs.LRLI(phylogeny, phylogeny.mat, all.matched.pairs, data.mat)
model3<-lm(ds$Log10.Rate~ds$Log10.Interval)

if (plot == "y")
{
	plot(ds$Log10.Rate~ds$Log10.Interval,ylab="Log10(Rate)",xlab="Log10(Interval)",cex.lab=1.3)
	abline(model3)
	par(mfrow=c(1,1))
}

to.return<-list(model1,model2,model3)

return(to.return)

}

