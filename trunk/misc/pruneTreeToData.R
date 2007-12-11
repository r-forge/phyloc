pruneTreeToData <-
function(x,tree,NA.omit=TRUE,...) {
    treeTaxa <- tree$tip.label
    traitTaxa <- NULL
	if (is.vector(x)) {
		if (NA.omit) traitTaxa <- names(na.omit(x[tree$tip.label]))
		else traitTaxa <- names(x[tree$tip.label])
	}
	else if (is.matrix(x) || is.data.frame(x)) {
		if (NA.omit) traitTaxa <- names(na.omit(x[tree$tip.label,]))
		else traitTaxa <- row.names(x[tree$tip.label,])
	}	
    trimTaxa <- setdiff(treeTaxa, traitTaxa)
    if (length(trimTaxa) > 0) 
        drop.tip(tree, trimTaxa, ...)
    else tree
}
