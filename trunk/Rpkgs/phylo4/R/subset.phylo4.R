################
# subset phylo4
################
#
setMethod("subset", "phylo4",
  function(x,tip.include=NULL,tips.exclude=NULL,node.subtree=NULL,...) {
    if (!is.null(tips.include)) {
      return(prune(x,x@tip.label[-tips]))
    }
    
    if (!is.null(tips.exclude)) {
      return(prune(x,tips))
    }
    
    if (!is.null(node.subtree)) {	
      #code from David Orme's clade.members function in CAIC
      if(length(node.subtree>1)) {
        warning(">1 node number supplied - only the first will be used")
        node.subtree <- node.subtree[1]
      }
      if(!(node.subtree %in% 1:nNodes(x))) {
        stop("Node number supplied not present in phylogeny")
      }
      while(length(intN <- x[x>phylo4::nTips(x)]) > 0){		
        minIntN <- min(intN)
        childOfMinIntN <- with(phy, edge[,2][which(edge[,1]==minIntN)])
        x <- c(x[x != minIntN], childOfMinIntN)
      }
      prune(x,x@tip.label[-unique(x)])
    }
    
    return(x)

})



