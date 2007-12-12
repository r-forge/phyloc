"find.cherries" <-
function(phylogeny){

#Finds the n cherries in a fully binary phylogeny
#Returns an n x 2 matrix with the tip names

ancestral.nodes<-which(table(phylogeny$edge[, 1][phylogeny$edge[, 2] > 0])==2)
ancestors<-names(ancestral.nodes)
pairs<-matrix(NA,nrow=length(ancestral.nodes),ncol=2)

for (i in 1:length(ancestors))
{
	sisters<-phylogeny$edge[,2][which(phylogeny$edge[,1]==ancestors[i])]
	pairs[i,1]<-phylogeny$tip.label[as.numeric(sisters[1])]
	pairs[i,2]<-phylogeny$tip.label[as.numeric(sisters[2])]
}
pairs
}

