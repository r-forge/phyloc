## as.phylo.R (2007-03-05)

##     Conversion Among Tree Objects

## Copyright 2005-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

old2new.phylo <- function(phy)
{
    mode(phy$edge) <- "numeric"
    phy$Nnode <- -min(phy$edge)
    n <- length(phy$tip.label)
    NODES <- phy$edge < 0
    phy$edge[NODES] <- n - phy$edge[NODES]
    phy
}

new2old.phylo <- function(phy)
{
    NTIP <- length(phy$tip.label)
    NODES <- phy$edge > NTIP
    phy$edge[NODES] <- NTIP - phy$edge[NODES]
    mode(phy$edge) <- "character"
    phy$Nnode <- NULL
    phy
}

as.phylo <- function (x, ...)
{
    if (class(x) == "phylo") return(x)
    UseMethod("as.phylo")
}

as.phylo.hclust <- function(x, ...)
{
    N <- dim(x$merge)[1]
    edge <- matrix(NA, 2*N, 2)
    edge.length <- numeric(2*N)
    ## `node' gives the number of the node for the i-th row of x$merge
    node <- numeric(N)
    node[N] <- N + 2
    cur.nod <- N + 3
    j <- 1
    for (i in N:1) {
        edge[j:(j + 1), 1] <- node[i]
        for (l in 1:2) {
            k <- j + l - 1
            if (x$merge[i, l] > 0) {
                edge[k, 2] <- node[x$merge[i, l]] <- cur.nod
                cur.nod <- cur.nod + 1
                edge.length[k] <- x$height[i] - x$height[x$merge[i, l]]
            } else {
                edge[k, 2] <- -x$merge[i, l]
                edge.length[k] <- x$height[i]
            }
        }
        j <- j + 2
    }
    if (is.null(x$labels))
      x$labels <- as.character(1:(N + 1))
    obj <- list(edge = edge, edge.length = edge.length,
                tip.label = x$labels, Nnode = N)
    class(obj) <- "phylo"
    reorder(obj)
}

as.phylo.phylog <- function(x, ...)
{
    tr <- read.tree(text = x$tre)
    n <- length(tr$tip.label)
    edge.length <- numeric(dim(tr$edge)[1])
    term  <- which(tr$edge[, 2] <= n)
    inte  <- which(tr$edge[, 2] > n)
    edge.length[term] <- x$leaves[tr$tip.label]
    edge.length[inte] <- x$nodes[tr$node.label][-1]
    tr$edge.length <- edge.length
    if (x$nodes["Root"] != 0) {
        tr$edge.root <- x$nodes["Root"]
        names(tr$edge.root) <- NULL
    }
    tr
}

as.hclust.phylo <- function(x, ...)
{
    if (!is.ultrametric(x)) stop("the tree is not ultrametric")
    if (!is.binary.tree(x)) stop("the tree is not binary")
    n <- length(x$tip.label)
    bt <- rev(branching.times(x))
    N <- length(bt)
    nm <- as.numeric(names(bt))
    merge <- matrix(NA, N, 2)
    for (i in 1:N) {
        ind <- which(x$edge[, 1] == nm[i])
        for (k in 1:2)
          merge[i, k] <- if (x$edge[ind[k], 2] <= n) -x$edge[ind[k], 2]
          else which(nm == x$edge[ind[k], 2])
    }
    names(bt) <- NULL
    obj <- list(merge = merge, height = bt, order = 1:(N + 1),
                labels = x$tip.label, call = match.call(),
                method = "unknown")
    class(obj) <- "hclust"
    obj
}
