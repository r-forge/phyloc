plot.multiPhylo <- function(x, layout = 1, ...)
{
    if (layout > 1)
      layout(matrix(1:layout, ceiling(sqrt(layout)), byrow = TRUE))
    if (!par("ask")) {
        par(ask = TRUE)
        on.exit(par(ask = FALSE))
    }
    ## assumes that this has only tip labels....
    LABELS <- attr(x, "Labels")
    for (i in x) {
        if (!is.null(LABELS)) i$tip.label <- LABELS
        plot(i, ...)
    }
}
