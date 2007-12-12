"root.to.tip" <-
function(phylogeny, tip.name){
#Compute a single root-to-tip value for taxon tip.name
n.label<-which(phylogeny$tip.label==tip.name)
n.edge<-which(phylogeny$edge[,2]==as.character(n.label))
sum.so.far<-phylogeny$changes[n.edge]
while(phylogeny$edge[n.edge,1]!="-1")
{
	n.edge<-which(phylogeny$edge[,2]==phylogeny$edge[n.edge,1])
	sum.so.far<-sum.so.far+phylogeny$changes[n.edge]
}
return(sum.so.far)
}

