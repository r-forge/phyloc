"matched.pairs.DI" <-
function(phylogeny, phylogeny.matrix, pairs, data.matrix){

#Computes differences in a trait and intervals for phylogenetically independent matched pairs
#'phylogeny' is the phylogeny
#'phylogeny.matrix' is a matrix of PD as constructed by cophenetic
#'pairs' is an n x 2 matrix of tip labels
#'data.matrix' is a matrix of data differences as constructed by as.matrix(dist())
#Returns a data frame with four columns, see code for details.


#10/5/06: Created from matched.pairs.LRLI, after bug fixes to that function

indices<-matrix(match(pairs[,], phylogeny$tip.label),nrow=length(pairs[,1]),ncol=2)
PD<-rep(NA,length(pairs[,1]))
X.diff<-rep(NA,length(pairs[,1]))

for (i in 1:length(PD))
{
	PD[i]<-phylogeny.matrix[indices[i,1],indices[i,2]]
	X.diff[i]<-data.matrix[indices[i,1],indices[i,2]]
}

X.diff[X.diff==0]<-min(X.diff[X.diff>0])/10 #Hack

to.return<-data.frame(pairs[,1],pairs[,2],X.diff,PD)
names(to.return)<-c("Tip.1","Tip.2","Difference","Interval")


return(to.return)
}

