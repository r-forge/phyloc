"get.root.to.tip" <-
function(phylogeny){
#Compute root to tip sum of character changes for each tip
r2t.changes<-rep(NA,length(phylogeny$tip.label))
for(i in 1:length(phylogeny$tip.label))
{	
	r2t.changes[i]<-root.to.tip(phylogeny,phylogeny$tip.label[i])
}
return(r2t.changes)
}

