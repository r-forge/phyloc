################
# subset phylo4
################
#

setGeneric("subset",signature(...))

setMethod("subset", "phylo4",
  function(object,tips.include=NULL,tips.exclude=NULL,node.subtree=NULL){

if (!is.null(tips.include)) {
  return(prune(object,object@tip.label[-tips]))
}

if (!is.null(tips.exclude)) {
  return(prune(object,tips))
}

if (!is.null(node.subtree)) {	
  #code from David Orme's clade.members function in CAIC
  if(length(node.subtree>1)) {
    warning(">1 node number supplied - only the first will be used")
    node.subtree <- node.subtree[1]
  }
  if(!(node.subtree %in% 1:nNodes(object))) {
    stop("Node number supplied not present in phylogeny")
  }
	while(length(intN <- object[object>phylo4::nTips(object)]) > 0){		
		minIntN <- min(intN)
		childOfMinIntN <- with(phy, edge[,2][which(edge[,1]==minIntN)])
		object <- c(object[object != minIntN], childOfMinIntN)
	}
  prune(object,object@tip.label[-unique(object)])
}

return(object)

})
