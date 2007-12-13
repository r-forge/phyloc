plot.phylo.coor <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, font = 3, cex = par("cex"), adj = NULL, srt = 0, 
    no.margin = FALSE, root.edge = FALSE, label.offset = 0, underscore = FALSE, 
    x.lim = NULL, y.lim = NULL, direction = "rightwards", lab4ut = "horizontal", 
    tip.color = "black", ...) 
{
    Ntip <- length(x$tip.label)
    if (Ntip == 1) 
        stop("found only one tip in the tree!")
    Nedge <- dim(x$edge)[1]
    if (any(tabulate(x$edge[, 1]) == 1)) 
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles().")
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
        "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
        "upwards", "downwards"))
    if (is.null(x$edge.length)) 
        use.edge.length <- FALSE
    if (type == "unrooted" || !use.edge.length) 
        root.edge <- FALSE
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    if (phyloORclado) {
        if (!is.null(attr(x, "order"))) 
            if (attr(x, "order") == "pruningwise") 
                x <- reorder(x)
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    edge.color <- rep(edge.color, length.out = Nedge)
    edge.width <- rep(edge.width, length.out = Nedge)
    xe <- x$edge
    x <- reorder(x, order = "pruningwise")
    ereorder <- match(x$edge[, 2], xe[, 2])
    edge.color <- edge.color[ereorder]
    edge.width <- edge.width[ereorder]
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) 
                node.pos <- 2
        }
        if (node.pos == 1) 
            yy <- .C("node_height", as.integer(Ntip), as.integer(Nnode), 
                as.integer(x$edge[, 1]), as.integer(x$edge[, 
                  2]), as.integer(Nedge), as.double(yy), DUP = FALSE, 
                PACKAGE = "ape")[[6]]
        else {
            ans <- .C("node_height_clado", as.integer(Ntip), 
                as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[, 
                  2]), as.integer(Nedge), double(Ntip + Nnode), 
                as.double(yy), DUP = FALSE, PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .C("node_depth", as.integer(Ntip), as.integer(Nnode), 
                  as.integer(x$edge[, 1]), as.integer(x$edge[, 
                    2]), as.integer(Nedge), double(Ntip + Nnode), 
                  DUP = FALSE, PACKAGE = "ape")[[6]] - 1
            xx <- max(xx) - xx
        }
        else {
            xx <- .C("node_depth_edgelength", as.integer(Ntip), 
                as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[, 
                  2]), as.integer(Nedge), as.double(x$edge.length), 
                double(Ntip + Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
        }
    }
    if (type == "fan") {
        TIPS <- xe[which(xe[, 2] <= Ntip), 2]
        xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
        theta <- double(Ntip)
        theta[TIPS] <- xx
        theta <- c(theta, numeric(Nnode))
        theta <- .C("node_height", as.integer(Ntip), as.integer(Nnode), 
            as.integer(x$edge[, 1]), as.integer(x$edge[, 2]), 
            as.integer(Nedge), theta, DUP = FALSE, PACKAGE = "ape")[[6]]
        if (use.edge.length) {
            r <- .C("node_depth_edgelength", as.integer(Ntip), 
                as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[, 
                  2]), as.integer(Nedge), as.double(x$edge.length), 
                double(Ntip + Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
        }
        else {
            r <- .C("node_depth", as.integer(Ntip), as.integer(Nnode), 
                as.integer(x$edge[, 1]), as.integer(x$edge[, 
                  2]), as.integer(Nedge), double(Ntip + Nnode), 
                DUP = FALSE, PACKAGE = "ape")[[6]]
            r <- 1/r
        }
        xx <- r * cos(theta)
        yy <- r * sin(theta)
    }
    if (type == "unrooted") {
        XY <- if (use.edge.length) 
            unrooted.xy(Ntip, Nnode, x$edge, x$edge.length)
        else unrooted.xy(Ntip, Nnode, x$edge, rep(1, Nedge))
        xx <- XY$M[, 1] - min(XY$M[, 1])
        yy <- XY$M[, 2] - min(XY$M[, 2])
    }
    if (type == "radial") {
        X <- .C("node_depth", as.integer(Ntip), as.integer(Nnode), 
            as.integer(x$edge[, 1]), as.integer(x$edge[, 2]), 
            as.integer(Nedge), double(Ntip + Nnode), DUP = FALSE, 
            PACKAGE = "ape")[[6]]
        X[X == 1] <- 0
        X <- 1 - X/Ntip
        yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
        Y <- .C("node_height", as.integer(Ntip), as.integer(Nnode), 
            as.integer(x$edge[, 1]), as.integer(x$edge[, 2]), 
            as.integer(Nedge), as.double(yy), DUP = FALSE, PACKAGE = "ape")[[6]]
        xx <- X * cos(Y)
        yy <- X * sin(Y)
    }
    if (phyloORclado && direction != "rightwards") {
        if (direction == "leftwards") {
            xx <- -xx
            xx <- xx - min(xx)
        }
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
            if (direction == "downwards") {
                yy <- -yy
                yy <- yy - min(yy)
            }
        }
    }
    if (phyloORclado && root.edge) {
        if (direction == "rightwards") 
            xx <- xx + x$root.edge
        if (direction == "upwards") 
            yy <- yy + x$root.edge
    }
    if (no.margin) 
        par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                tmp <- if (show.tip.label) 
                  nchar(x$tip.label) * 0.018 * max(xx) * cex
                else 0
                x.lim[2] <- if (direction == "leftwards") 
                  max(xx[ROOT] + tmp)
                else max(xx[1:Ntip] + tmp)
            }
            else x.lim <- c(1, Ntip)
        }
        if (type == "fan") {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
            }
            else x.lim <- c(min(xx), max(xx))
        }
        if (type == "unrooted") {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            }
            else x.lim <- c(0, max(xx))
        }
        if (type == "radial") {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            }
            else x.lim <- c(-1, 1)
        }
    }
    else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            x.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.03 * cex)
            else -1
    }
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
                y.lim <- c(1, Ntip)
            else {
                y.lim <- c(0, NA)
                tmp <- if (show.tip.label) 
                  nchar(x$tip.label) * 0.018 * max(yy) * cex
                else 0
                y.lim[2] <- if (direction == "downwards") 
                  max(yy[ROOT] + tmp)
                else max(yy[1:Ntip] + tmp)
            }
        }
        if (type == "fan") {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            }
            else y.lim <- c(min(yy), max(yy))
        }
        if (type == "unrooted") {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            }
            else y.lim <- c(0, max(yy))
        }
        if (type == "radial") {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            }
            else y.lim <- c(-1, 1)
        }
    }
    else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            y.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
            else -1
    }
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
            x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
            y.lim[2] <- y.lim[2] + x$root.edge
    }
    #########################################HERE
       if (is.null(adj)) 
        adj <- if (phyloORclado && direction == "leftwards") 
            1
        else 0
    if (phyloORclado) {
        MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
        if (direction == "rightwards") {
            lox <- label.offset + MAXSTRING * 1.05 * adj
            loy <- 0
        }
        if (direction == "leftwards") {
            lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
            loy <- 0
            xx <- xx + MAXSTRING
        }
        if (!horizontal) {
            psr <- par("usr")
            MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                psr[1])
            loy <- label.offset + MAXSTRING * 1.05 * adj
            lox <- 0
            srt <- 90 + srt
            if (direction == "downwards") {
                loy <- -loy
                yy <- yy + MAXSTRING
                srt <- 180 + srt
            }
        }
    }
      
    
    if (root.edge) 
        switch(direction, rightwards = segments(0, yy[ROOT], 
            x$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
            yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]), upwards = segments(xx[ROOT], 
            0, xx[ROOT], x$root.edge), downwards = segments(xx[ROOT], 
            yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge))
    if (show.tip.label) {
        if (!underscore) 
            x$tip.label <- gsub("_", " ", x$tip.label)
               if (type == "unrooted") {
            if (lab4ut == "horizontal") {
                y.adj <- x.adj <- numeric(Ntip)
                sel <- abs(XY$axe) > 0.75 * pi
                x.adj[sel] <- -strwidth(x$tip.label)[sel] * 1.05
                sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
                  pi
                x.adj[sel] <- -strwidth(x$tip.label)[sel] * (2 * 
                  abs(XY$axe)[sel]/pi - 0.5)
                sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
                y.adj[sel] <- strheight(x$tip.label)[sel]/2
                sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
                y.adj[sel] <- -strheight(x$tip.label)[sel] * 
                  0.75
                
            }
            else {
                adj <- as.numeric(abs(XY$axe) > pi/2)
                srt <- 180 * XY$axe/pi
                srt[as.logical(adj)] <- srt[as.logical(adj)] - 
                  180
                sel <- srt > -1e-06 & srt < 0
                if (any(sel)) 
                  srt[sel] <- 0
                
            }
        }
        if (type %in% c("fan", "radial")) {
            xx.scaled <- xx[1:Ntip]
            if (type == "fan") {
                maxx <- max(xx.scaled)
                if (maxx > 1) 
                  xx.scaled <- xx.scaled/maxx
            }
            angle <- acos(xx.scaled) * 180/pi
            s1 <- angle > 90 & yy[1:Ntip] > 0
            s2 <- angle < 90 & yy[1:Ntip] < 0
            s3 <- angle > 90 & yy[1:Ntip] < 0
            angle[s1] <- angle[s1] + 180
            angle[s2] <- -angle[s2]
            angle[s3] <- 180 - angle[s3]
            adj <- numeric(Ntip)
            adj[xx[1:Ntip] < 0] <- 1
           
        }
    }
   
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
    .last_plot.phylo <<- c(L, list(edge = xe, xx = xx, yy = yy))
    invisible(L)
    ##lines for developement
    coor<-matrix(ncol=2, nrow=Ntip+Nnode)
    coor[,1]<-xx
    coor[,2]<-yy
    #row.names(coor)<-x$tip.label
    return(coor)
    #return(x.lim)
}
