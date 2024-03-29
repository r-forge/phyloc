## ltt.plot.R (2007-12-11)

##    Lineages Through Time Plot

## Copyright 2002-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ltt.plot <- function(phy, xlab = "Time", ylab = "N", ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    time <- sort(branching.times(phy))
    N <- 1:(length(time) + 1)
    plot(-c(rev(time), 0), N, xlab = xlab, ylab = ylab,
         xaxs = "r", yaxs = "r", type = "S", ...)
}

ltt.lines <- function(phy, ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    time <- sort(branching.times(phy))
    N <- 1:(length(time) + 1)
    lines(-c(rev(time), 0), N, type = "S", ...)
}

mltt.plot <- function(phy, ..., dcol = TRUE, dlty = FALSE, legend = TRUE,
                      xlab = "Time", ylab = "N")
{
    if (!class(phy) %in% c("phylo", "multiPhylo"))
        stop("tree(s) are not of appropriate class")

    ltt.xy <- function(phy) {
        x <- -c(rev(sort(branching.times(phy))), 0)
        names(x) <- NULL
        y <- 1:length(x)
        cbind(x, y)
    }
    if (class(phy) == "phylo") {
        TREES <- list(ltt.xy(phy))
        names(TREES) <- deparse(substitute(phy))
    } else { # class(phy) == "multiPhylo"
        TREES <<- TREES
        for (i in 1:length(TREES))
            phy[[i]]$tip.label <- attr(phy, "Labels")
        TREES <- lapply(phy, ltt.xy)
#        TREES <<- TREES
        names(TREES) <- names(phy)
    }
    dts <- list(...)
    if (length(dts)) {
        mc <- as.character(match.call())[-(1:2)]
        nms <- mc[1:length(dts)]
        for (i in 1:length(dts)) {
            if (length(class(dts[[i]])) == 1) {
                a <- list(ltt.xy(dts[[i]]))
                names(a) <- nms[i]
            } else {
                a <- lapply(dts[[i]], ltt.xy)
                names(a) <- names(dts[[i]])
            }
            TREES <- c(TREES, a)
        }
    }
    n <- length(TREES)
    xl <- c(min(unlist(lapply(TREES, function(x) min(x[, 1])))), 0)
    yl <- c(1, max(unlist(lapply(TREES, function(x) max(x[, 2])))))

    plot(0, 0, type = "n", xlim = xl, ylim = yl, xaxs = "r", yaxs = "r",
         xlab = xlab, ylab = ylab)

    lty <- if (!dlty) rep(1, n) else 1:n
    col <- if (!dcol) rep(1, n) else topo.colors(n)

    for (i in 1:n)
      lines(TREES[[i]], col = col[i], lty = lty[i], type = "S")

    if (legend)
      legend(xl[1], yl[2], legend = names(TREES),
             lty = lty, col = col, bty = "n")
}
