
## REQUIRED for all trees
check_edgestruc <- function(object) {
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
                       missing.tip=c("stop","OK","warn"),
                       missing.node=c("OK","warn","stop"),
                       missing.edge=c("OK","warn","stop"),
                       extra.tip=c("stop","OK","warn"),
                       extra.node=c("OK","warn","stop"),
                       extra.edge=c("OK","warn","stop")) {
  ## do we also have a set of these for names?
  
  ## tip default: must match exactly (missing=="stop", extra=="stop")
  ## node default:

  ## this sets the default
  missing.tip <- match.arg(missing.tip)

}
                       
                       

