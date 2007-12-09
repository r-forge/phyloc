`evolve.brownian` <-
function(phy,value=0,var=1) {
	x = as.matrix(evolve.phylo(phy,value,var)$tip.character)
	row.names(x) = phy$tip.label
	return(x)
	}

