"get.branch.length" <-
function(phylogeny, tip.name){

#Returns the edge length leading to a particular tip in the phylogeny

edge.label<-which(phylogeny$tip.label==tip.name)
edge.index<-which(phylogeny$edge[,2]==as.character(edge.label))

return(phylogeny$edge.length[edge.index])
}

