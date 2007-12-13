"clade.matrix" <-
function(phy){

	# returns the phylogeny as a table showing clade
	# membership below each node. 

	# check the object is a phylogeny
	if(class(phy) != "phylo") stop("Phylogeny required")
	
	# get the number of tips and nodes
	nb.tips <- max(as.numeric(phy$edge))
	nb.nodes <- -min(as.numeric(phy$edge))
	
	# set-up the clade matrix
	mat.names <- list(edges= c(1:nb.tips, (-nb.nodes):-1), tips=1:nb.tips)
	clade.mat <- matrix(0, nrow=nb.tips+nb.nodes, ncol=nb.tips, dimnames=mat.names)
	
	# the diagonal from [1,1] are the tips
	diag(clade.mat) <- 1
	
	# now deal with internals
	node.members <- all.clades(phy)

	node.id <- names(node.members)
	for(rows in 1:length(node.id)){
		clade.mat[node.id[rows], node.members[rows][[1]]] <- 1}
	
	# get edge lengths into correct order, inserting root edge
	if(is.null(phy$root.edge)){
		edge.len <- c(phy$edge.length,0)
		} else {
		edge.len <- c(phy$edge.length,phy$root.edge)}
	names(edge.len) <- c(phy$edge[,2], "-1")
	edge.len <- edge.len[as.character(mat.names$edges)]
	
	RET <- list(clade.matrix=clade.mat, edge.length=edge.len, tip.label=phy$tip.label)
	class(RET) <- "clade.matrix"
	
	return(RET)
}

