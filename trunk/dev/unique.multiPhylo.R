unique.multiPhylo <- function(x, incomparables = FALSE,
                              use.edge.length = FALSE,
                              use.tip.label = TRUE, ...)
{
    n <- length(x)
    for (i in 1:n) x[[i]]$tip.label <- attr(x, "Labels")
    keep <- !logical(n)
    for (i in 2:n) {
        j <- 1
        while (j < i) {
            if (all.equal(x[[j]], x[[i]],
                          use.edge.length = use.edge.length,
                          use.tip.label = use.tip.label)) {
                keep[i] <- FALSE
                break
            }
            j <- j + 1
        }
    }
    structure(x[keep], class = "multiPhylo")
}
