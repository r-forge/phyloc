"clade.matrix" <-
function(phy){

    # OLD2NEW: CONVERTED
    
	# returns the phylogeny as a table showing clade
	# membership below each node. 

	# check the object is a phylogeny
	if(class(phy) != "phylo") stop("Phylogeny required")
	
	# get the number of tips and nodes
	nb.tips <- max(phy$edge) - phy$Nnode
	nb.nodes <- phy$Nnode
	nb.all <- nb.tips + nb.nodes
	
	# set-up the clade matrix
	mat.names <- list(edges= c(1:nb.all), tips=1:nb.tips)
	clade.mat <- matrix(0, nrow=nb.all, ncol=nb.tips, dimnames=mat.names)
	
	# the diagonal from [1,1] are the tips
	diag(clade.mat) <- 1
	
	# now deal with internals
	node.members <- all.clades(phy)

	node.id <- names(node.members)
	for(rows in seq(along=node.id)){
		clade.mat[node.id[rows], node.members[rows][[1]]] <- 1}
	
	# get edge lengths into correct order, inserting root edge
	if(is.null(phy$root.edge)){
		edge.len <- c(phy$edge.length,0)
		} else {
		edge.len <- c(phy$edge.length,phy$root.edge)}
	
	names(edge.len) <- c(phy$edge[,2], nb.tips + 1)
	edge.len <- edge.len[as.character(mat.names$edges)]
	
	RET <- list(clade.matrix=clade.mat, edge.length=edge.len, tip.label=phy$tip.label)
	class(RET) <- "clade.matrix"
	
	return(RET)
}

