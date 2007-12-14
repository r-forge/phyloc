################
# subset phylo4
################
#
setMethod("subset", "phylo4",
  function(x,tips.include=NULL,tips.exclude=NULL,node.subtree=NULL,...) {
    if (!is.null(tips.include)) {
      if (is.numeric(tips.include)) {
        tips.include <- x@tip.label[tips.include]
      }
      return(prune(x,x@tip.label[-na.omit(match(x@tip.label,tips.include))]))
    }
    
    if (!is.null(tips.exclude)) {
      return(prune(x,tips.exclude))
    }
    
    if (!is.null(node.subtree)) {
      #code from David Orme's clade.members function in CAIC
      if(length(node.subtree>1)) {
        warning(">1 node number supplied - only the first will be used")
        node.subtree <- node.subtree[1]
      }
      if(!(node.subtree %in% 1:(phylo4::nTips(x)+nNodes(x)))) {
        stop("Node number supplied not present in phylogeny")
      }
    }
    
    return(x)
  })



setMethod("subset", "phylo", function(x,tips.include=NULL,tips.exclude=NULL,node.subtree=NULL,...) {

  x <- as(x,"phylo")
  res <- subset(x,tips.include=NULL,tips.exclude=NULL,node.subtree=NULL,...)
  return(res)
  })



#
# PLEASE CODE ME 
#
setMethod("subset", "phylo4d", function(x,tips.include=NULL,tips.exclude=NULL,node.subtree=NULL,...) {

# place commercial here.
})

