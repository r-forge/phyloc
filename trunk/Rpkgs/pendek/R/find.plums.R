"find.plums" <-
function(phylogeny){

#Finds the plums in a phylogeny that can contain polytomies
#A 'plum' (a neologism) is a clade composed of the descendants of a node, all of whose
#descendants are tips (i.e., a bit like cherries, but can be bigger).
#Returns a list of vectors, each vector holding the tip names from one plum

all.nodes<-unique(phylogeny$edge[,1])
nodes.with.nodes.descending<-names(table(phylogeny$edge[, 1][phylogeny$edge[, 2] < 0]))
plum.nodes<-setdiff(all.nodes,nodes.with.nodes.descending)

plum.list<-NULL
for (i in 1:length(plum.nodes))
{
	members<-phylogeny$tip.label[as.numeric(phylogeny$edge[,2][which(phylogeny$edge[,1]==plum.nodes[i])])]
	plum<-list(node=plum.nodes[i],tips=members)
	plum.list<-c(plum.list,list(plum))
}
plum.list
}

