'pic.NA' <-
function(x,tree,...) {
	
  nNodes <- dim(tree$edge)[1] + 1
  nIntNodes <- nNodes - length(tree$tip.label)
	tree$node.label <- paste("node",1:nIntNodes,sep="")
  output <- data.frame(rep(NA,nIntNodes),row.names=tree$node.label)
  for (i in 1:dim(x)[2]) {
    trait <- x[,i]
    names(trait) <- row.names(x)
    pruned <- pruneMissing(trait,tree)
    result <- pic(pruned$data,pruned$tree,...)
    names(result) <- pruned$tree$node.label
    output[,i] <- result[match(row.names(output),names(result))]
  }
  if (!is.null(colnames(x))) colnames(output) <- colnames(x)
  output

}