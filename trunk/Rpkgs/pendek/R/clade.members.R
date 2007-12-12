"clade.members" <-
function(x, phyl, tip.labels=FALSE){
	
	# returns a vector of the tips that subtend from an
	# identified node
	if(class(phyl) != "phylo") stop("Phylogeny required")
	
	if(!(x %in% as.numeric(phyl$edge))) stop("Node not in range for phylogeny")
	
	while(any(x < 0)){
		minx.loc <- which(x == min(x))
		minx <- min(x)
		
		x <- x[-minx.loc]
		newx <- as.numeric(phyl$edge[,2][which(phyl$edge[,1]==minx)])
		x <- c(x, newx)
	}
	
	x <- unique(x)
	
	if(tip.labels) {
		return(phyl$tip.label[x])
	} else {
		return(x)
	}
}

