print.multi.tree<-function (x, details = FALSE, printlen=6,...) 
{
    N <- length(x)
    nb.tip <- Ntip (x[[1]])
    nb.node<- Nnode(x[[1]])

   cat(paste("\n",N,"phylogenetic trees with", nb.tip, "tips and", 
        nb.node, "internal nodes.\n\n"))
    cat("Tip labels:\n")
    if (nb.tip > printlen) {
        cat(paste("\t", paste(x[[1]]$tip.label[1:printlen], collapse = ", "), 
            ", ...\n", sep = ""))
    }
    else print(x[[1]]$tip.label)
    if (!is.null(x[[1]]$node.label)) {
        cat("\tNode labels:\n")
        if (nb.node > printlen) {
            cat(paste("\t", paste(x[[1]]$node.label[1:printlen], collapse = ", "), 
                ",...\n", sep = ""))
        }
        else print(x[[1]]$node.label)
    }
    rlab <- if (is.rooted(x$tree1)) 
        "Rooted"
    else "Unrooted"
    cat("\n", rlab, "; ", sep = "")
    blen <- if (is.null(x[[1]]$edge.length)) 
        "no branch lengths."
    else "includes branch lengths."
    cat(blen, "\n", sep = "")



    if (details) 
        for (i in 1:N) cat("tree", i, ":", length(x[[i]]$tip.label), 
            "tips\n")
    cat("\n")
}
