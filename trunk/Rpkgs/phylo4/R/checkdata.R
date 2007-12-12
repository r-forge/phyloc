
## REQUIRED for all trees
check_phylo4 <- function(object) {
	## check consistency of edges etc.
	N <- nrow(object@edge)
	if (hasEdgeLength(object) && length(object@edge.length) != N)
		return("edge lengths do not match number of edges")
	if (length(object@tip.label)+object@Nnode-1 != N)
		return("number of tip labels not consistent with number of edges and nodes")
	return(TRUE)
}


## VERY limited: required of all phylo4d objects
check_data0 <- function(object) {
	## check consistency 
	return(TRUE)
}


check_data <- function(object,
		use.tip.names=c(TRUE,FALSE),
		missing.tip.data=c("stop","OK","warn"),
		extra.tip.data=c("stop","OK","warn"),
		missing.tip.names=c("stop","OK","warn"),
		use.node.names=c(FALSE,TRUE),
		missing.node.data=c("OK","warn","stop"),		
		extra.node.data=c("OK","warn","stop"),												
		missing.node.names=c("OK","warn","stop"),
		use.edge.names=c(FALSE,TRUE),
		missing.edge.data=c("OK","warn","stop"),
		extra.edge.data=c("OK","warn","stop"),													 
		missing.edge.names=c("OK","warn","stop"))											 
{
		
	## tip default: use names, require names, must match exactly
	##(missing.data=="stop", extra.data=="stop", missing.names=="stop")
	use.tip.names=match.arg(use.tip.names)
	missing.tip.data <- match.arg(missing.tip)
	extra.tip.data <- match.arg(extra.tip)
	missing.tip.names <- match.arg(missing.tip.names)
	
	## node default: don't use node names, don't require names, do not need to match exactly
	use.node.names=match.arg(use.node.names)
	missing.node.data <- match.arg(missing.node)
	extra.node.data <- match.arg(extra.node)
	missing.node.names <- match.arg(missing.tip.names)
	
	## edge default: don't use edge names, don't require names, do not need to match exactly
	use.node.names=match.arg(use.edge.names)
	missing.edge.data <- match.arg(missing.edge)
	extra.edge.data <- match.arg(extra.edge)
	missing.edge.names <- match.arg(missing.tip.names)
	
	## clean empty names - if data are imported as 'all' but only tips or nodes have names, clean up
	if (all(row.names(object@tipdata)==""))
		row.names(object@tipdata) <- NULL
	if (all(row.names(object@nodedata)==""))
		row.names(object@nodedata) <- NULL
	
	## for each set of data, check for names, missing and extra data and take appropriate actions
	
	## tip data checks
	## if tipdata exist
	if (!is.null(object@tipdata)) {
		## if we want to use tip.names
		if (use.tip.names) {
			## check for names
			if (!is.null(row.names(object@tipdata))) {
				## check for missing or extra tip data (relative to tree taxa)
				if (setequal(row.names(object@tipdata), object@tip.label)) {
					##names are perfect match - ok
					return(TRUE)
				}
				else {
					#we know the tree taxa and tipdata taxa are not a perfect match
					#if tree taxa are subset of tipdata, check missing.tip arg and act accordingly
					if (all(row.names(object@tipdata) %in% object@tip.label)) {
						#we know it's not an exact match - we have missing.tip.data - take action
						#stop
						if (missing.tip.data == "stop") {
							return("Tip data names do not exactly match phylo4 tip labels.")
						}
						#warn
						else if (missing.tip.data == "warn") {
							warning("Tip data names do not exactly match phylo4 labels.")
							return(TRUE)
						}
						#else ok
					}
					#if tipdata taxa are subset of tree taxa, check extra.tip arg and act accordingly
					if (all(object@tip.label %in% row.names(object@tipdata))) {
						#we know it's not an exact match - we have extra.tip.data - take action
						#stop
						if (missing.tip.data == "stop") {
							return("Tip data are a superset of phylo4 tips.")
						}
						else if (missing.tip.data == "warn") {
							warning("Tip data are a superset of phylo4 tips.")
							return(TRUE)
						}
						#else ok
					}
					return(TRUE)
				}
			}
			else {
				#no tip.names
				if (missing.tip.names == "stop") {
					return("Tip data do not have names.")
				}
				else if (missing.tip.names == "warn") {
					warning("Tip data do not have names.")
					return(TRUE)
				}
				#don't use tip names or attempt to sort - but check to make sure dimensions match
				if (!(nTips(tree)==length(object@tipdata))) {
					return("Tip data do not have names and do not match number of phylo4 tips.")
				}
			}
		}
	}
	else
	{
		#don't use tip names or attempt to sort - but check to make sure dimensions match
		if (!(nTips(tree)==length(object@tipdata))) {
			return("Tip data do not have names and do not match number of phylo4 tips.")
		}
	}

## node data checks
	## if tipdata exist
	if (!is.null(object@nodedata)) {
		## if we want to use node.names
		if (use.node.names) {
			## check for names
			if (!is.null(row.names(object@nodedata))) {
				## check for missing or extra node data (relative to tree taxa)
				if (setequal(row.names(object@nodedata), object@node.label)) {
					##names are perfect match - ok
					return(TRUE)
				}
				else {
					#we know the tree taxa and nodedata taxa are not a perfect match
					#if tree taxa are subset of nodedata, check missing.node arg and act accordingly
					if (all(row.names(object@nodedata) %in% object@node.label)) {
						#we know it's not an exact match - we have missing.node.data - take action
						#stop
						if (missing.node.data == "stop") {
							return("Node data names do not exactly match phylo4 node labels.")
						}
						#warn
						else if (missing.node.data == "warn") {
							warning("Node data names do not exactly match phylo4 labels.")
							return(TRUE)
						}
						#else ok
					}
					#if nodedata taxa are subset of tree taxa, check extra.node arg and act accordingly
					if (all(object@node.label %in% row.names(object@nodedata))) {
						#we know it's not an exact match - we have extra.node.data - take action
						#stop
						if (missing.node.data == "stop") {
							return("Node data are a superset of phylo4 nodes.")
						}
						else if (missing.node.data == "warn") {
							warning("Node data are a superset of phylo4 nodes.")
							return(TRUE)
						}
						#else ok
					}
					return(TRUE)
				}
			}
			else {
				#no node.names
				if (missing.node.names == "stop") {
					return("Node data do not have names.")
				}
				else if (missing.node.names == "warn") {
					warning("Node data do not have names.")
					return(TRUE)
				}
				#don't use node names or attempt to sort - but check to make sure dimensions match
				if (!(nNodes(tree)==length(object@nodedata))) {
					return("Node data do not have names and do not match number of phylo4 nodes.")
				}
			}
		}
	}
	else
	{
		#don't use node names or attempt to sort - but check to make sure dimensions match
		if (!(nNodes(tree)==length(object@nodedata))) {
			return("Tip data do not have names and do not match number of phylo4 nodes.")
		}
	}

## edge data checks
## TODO check this section - do edge.label and such exist? do we just use node.label? or ignore?

}


attach_data <- function(object,
		use.tip.names=c(TRUE,FALSE),
		use.node.names=c(FALSE,TRUE),
		use.edge.names=c(FALSE,TRUE))							 
{
	
	## assumes data have already been checked by check_data
	use.tip.names=match.arg(use.tip.names)
	use.node.names=match.arg(use.node.names)
	use.node.names=match.arg(use.edge.names)
	
	## clean empty names - if data are imported as 'all' but only tips or nodes have names, clean up
	if (all(row.names(object@tipdata)==""))
		row.names(object@tipdata) <- NULL
	if (all(row.names(object@nodedata)==""))
		row.names(object@nodedata) <- NULL
	
	## for each set of data, take appropriate actions
	
	## tip data operations:
	## if tipdata exist
	if (!is.null(object@tipdata)) {
		## if we want to use tip.names
		if (use.tip.names) {
			tipdata <- tipdata[match(row.names(object@tipdata),object@tip.label),]
			row.names(tipdata) <- object@tip.label
		}
	}
	
	## node data operations
	if (!is.null(object@nodedata)) {
		## if we want to use tip.names
		if (use.tip.names) {
			nodedata <- nodedata[match(row.names(object@nodedata),object@node.label),]
			row.names(nodedata) <- object@node.label
		}
	}

	## edge data operations
	if (!is.null(object@edgedata)) {
		## if we want to use tip.names
		if (use.edge.names) {
			edgedata <- edgedata[match(row.names(object@edgedata),object@edge.label),]
			row.names(tipdata) <- object@tip.label
		}
	}
	
	return(object)
	
}
