"all.clades" <-
function(phyl, tips=FALSE, tip.labels=FALSE){

	# returns a list of vectors showing the tips
	# subtending from each node in the tree
	if(class(phyl) != "phylo") stop("Phylogeny required")

	nodes <- sort(unique(as.numeric(phyl$edge)))
	
	if(!tips) nodes <- nodes[nodes < 0]
	
	clade.list <- mapply(clade.members, nodes, MoreArgs=list(phyl=phyl, tip.labels=tip.labels))
	names(clade.list) <- nodes
	
	return(clade.list)
}

