"matched.pairs.LRLI" <-
function(phylogeny, phylogeny.matrix, pairs, data.matrix, plot="n"){

#Plots and models LRLI based on phylogenetically independent matched pairs
#'phylogeny' is the phylogeny
#'phylogeny.matrix' is a matrix of PD as constructed by dist.phylo
#'pairs' is an n x 2 matric of tip labels
#'data.matrix' is a matrix of data differences as constructed by as.matrix(dist())
#Returns a data frame (used to be a plot and a model).

#10/5/06: Bug fix to prevent a rate of zero being found (which, when logged, causes error).
#See also all.LRLI

#10/5/06: Bug fix - was not previously passed the phylogeny, which it needed for the tip.labels

#10/5/06: Structure changed to return data frame rather than model.

indices<-matrix(match(pairs[,], phylogeny$tip.label),nrow=length(pairs[,1]),ncol=2)
PD<-rep(NA,length(pairs[,1]))
X.diff<-rep(NA,length(pairs[,1]))

for (i in 1:length(PD))
{
	PD[i]<-phylogeny.matrix[indices[i,1],indices[i,2]]
	X.diff[i]<-data.matrix[indices[i,1],indices[i,2]]
}

X.diff[X.diff==0]<-min(X.diff[X.diff>0])/10 #Hack
Log10.Rate<-log10(X.diff/PD)
Log10.Interval<-log10(PD)

to.return<-data.frame(pairs[,1],pairs[,2],Log10.Rate,Log10.Interval)
names(to.return)[1:2]<-c("Tip.1","Tip.2")

return(to.return)
}

