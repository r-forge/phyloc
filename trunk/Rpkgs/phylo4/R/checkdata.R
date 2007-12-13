
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


check_data <- function(object,
		use.tip.names=TRUE,
		missing.tip.data=c("fail","OK","warn"),
		extra.tip.data=c("fail","OK","warn"),
		missing.tip.names=c("fail","OK","warn"),
		use.node.names=FALSE,
		missing.node.data=c("OK","warn","fail"),		
		extra.node.data=c("OK","warn","fail"),												
		missing.node.names=c("OK","warn","fail"),...)							 

{
		
	## tip default: use names, require names, must match exactly
	#use.tip.names=match.arg(use.tip.names)
	missing.tip.data <- match.arg(missing.tip.data)
	extra.tip.data <- match.arg(extra.tip.data)
	missing.tip.names <- match.arg(missing.tip.names)
	
	## node default: don't use node names, don't require names, do not need to match exactly
	#use.node.names=match.arg(use.node.names)
	missing.node.data <- match.arg(missing.node.data)
	extra.node.data <- match.arg(extra.node.data)
	missing.node.names <- match.arg(missing.tip.names)
	
	## clean empty names - if data are imported as 'all' but only tips or nodes have names, clean up
	if (all(row.names(object@tip.data)==""))
		row.names(object@tip.data) <- NULL
	if (all(row.names(object@node.data)==""))
		row.names(object@node.data) <- NULL
	
	## for each set of data, check for names, missing and extra data and take appropriate actions
	
	## tip data checks
	## if tip.data exist
	if (!all(dim(object@tip.data)==0)) {
		## if we want to use tip.names
		if (use.tip.names) {
			## check for names
			if (!is.null(row.names(object@tip.data))) {
				## check for missing or extra tip data (relative to tree taxa)
				if (setequal(row.names(object@tip.data), object@tip.label)) {
					##names are perfect match - ok
					return(TRUE)
				}
				else {
					#we know the tree taxa and tip.data taxa are not a perfect match
					#if tree taxa are subset of tip.data, check missing.tip arg and act accordingly
					if (all(row.names(object@tip.data) %in% object@tip.label)) {
						#we know it's not an exact match - we have missing.tip.data - take action
						#fail
						if (missing.tip.data == "fail") {
							return("Tip data names do not exactly match phylo4 tip labels.")
						}
						#warn
						else if (missing.tip.data == "warn") {
							warning("Tip data names do not exactly match phylo4 labels.")
							return(TRUE)
						}
						#else ok
					}
					#if tip.data taxa are subset of tree taxa, check extra.tip arg and act accordingly
					if (all(object@tip.label %in% row.names(object@tip.data))) {
						#we know it's not an exact match - we have extra.tip.data - take action
						#fail
						if (missing.tip.data == "fail") {
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
				if (missing.tip.names == "fail") {
					return("Tip data do not have names.")
				}
				else if (missing.tip.names == "warn") {
					warning("Tip data do not have names.")
					return(TRUE)
				}
				#don't use tip names or attempt to sort - but check to make sure dimensions match
				if (!(phylo4::nTips(object)==length(object@tip.data))) {
					return("Tip data do not have names and do not match number of phylo4 tips.")
				}
			}
		}
	}
	else
	{
		#don't use tip names or attempt to sort - but check to make sure dimensions match
		if (!(phylo4::nTips(object)==length(object@tip.data))) {
			return("Tip data do not have names and do not match number of phylo4 tips.")
		}
	}

	## node data checks
	## if tip.data exist
	if (!all(dim(object@node.data)==0)) {
		## if we want to use node.names
		if (use.node.names) {
			## check for names
			if (!is.null(row.names(object@node.data))) {
				## check for missing or extra node data (relative to tree taxa)
				if (setequal(row.names(object@node.data), object@node.label)) {
					##names are perfect match - ok
					return(TRUE)
				}
				else {
					#we know the tree taxa and node.data taxa are not a perfect match
					#if tree taxa are subset of node.data, check missing.node arg and act accordingly
					if (all(row.names(object@node.data) %in% object@node.label)) {
						#we know it's not an exact match - we have missing.node.data - take action
						#fail
						if (missing.node.data == "fail") {
							stop("Node data names do not exactly match phylo4 node labels.")
						}
						#warn
						else if (missing.node.data == "warn") {
							warning("Node data names do not exactly match phylo4 labels.")
							return(TRUE)
						}
						#else ok
					}
					#if node.data taxa are subset of tree taxa, check extra.node arg and act accordingly
					if (all(object@node.label %in% row.names(object@node.data))) {
						#we know it's not an exact match - we have extra.node.data - take action
						#fail
						if (missing.node.data == "fail") {
							stop("Node data are a superset of phylo4 nodes.")
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
				if (missing.node.names == "fail") {
					return("Node data do not have names.")
				}
				else if (missing.node.names == "warn") {
					warning("Node data do not have names.")
					return(TRUE)
				}
				#don't use node names or attempt to sort - but check to make sure dimensions match
				if (!(nNodes(object)==length(object@node.data))) {
					stop("Node data do not have names and do not match number of phylo4 nodes.")
				}
			}
		}
	}
	else
	{
		#don't use node names or attempt to sort - but check to make sure dimensions match
		if (!(nNodes(object)==length(object@node.data))) {
			stop("Node data do not have names and do not match number of phylo4 nodes.")
		}
	}
}


attach_data <- function(object,
		use.tip.names=TRUE,
		use.node.names=FALSE,
		use.edge.names=FALSE,
		...)							 
{
	
	## assumes data have already been checked by check_data!
	
	## clean empty names - if data are imported as 'all' but only tips or nodes have names, clean up
	if (all(row.names(object@tip.data)==""))
		row.names(object@tip.data) <- NULL
	if (all(row.names(object@node.data)==""))
		row.names(object@node.data) <- NULL
	
	## for each set of data, take appropriate actions
	
	## tip data operations:
	## if tip.data exist
	if (!all(dim(object@tip.data)==0)) {
		## if we want to use tip.names
		if (use.tip.names) {
			object@tip.data <- object@tip.data[match(row.names(object@tip.data),object@tip.label),]
			row.names(object@tip.data) <- object@tip.label
		}
	}
	
	## node data operations
	if (!all(dim(object@node.data)==0)) {
		## if we want to use tip.names
		if (use.node.names) {
			object@node.data <- object@node.data[match(row.names(object@node.data),object@node.label),]
			row.names(object@node.data) <- object@node.label
		}
	}
	
	return(object)
	
}
