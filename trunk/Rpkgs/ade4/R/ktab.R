########### is.ktab ###########
"is.ktab" <- function (x)
    inherits(x, "ktab")

########### [.ktab" ########### 
"[.ktab" <- function (x, selection) {
    blocks <- x$blo
    nblo <- length(blocks)
    if (is.logical(selection)) 
        selection <- which(selection)
    if (any(selection > nblo)) 
        stop("Non convenient selection")
    indica <- as.factor(rep(1:nblo, blocks))
    res <- unclass(x)[selection]
    cw <- x$cw
    cw <- split(cw, indica)
    cw <- unlist(cw[selection])
    res$cw <- cw
    res$lw <- x$lw
    blocks <- unlist(lapply(res, function(x) ncol(x)))
    nblo <- length(blocks)
    res$blo <- blocks
    ktab.util.addfactor(res) <- list(blocks, length(res$lw))
    res$call <- match.call()
    class(res) <- "ktab"
    return(res)
}


########### print.ktab ########### 
"print.ktab" <- function (x, ...) {
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    cat("class:", class(x), "\n")
    ntab <- length(x$blo)
    cat("\ntab number:  ", ntab, "\n")
    sumry <- array("", c(ntab, 3), list(1:ntab, c("data.frame", 
        "nrow", "ncol")))
    for (i in 1:ntab) {
        sumry[i, ] <- c(names(x)[i], nrow(x[[i]]), ncol(x[[i]]))
    }
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(4, 4), list((ntab + 1):(ntab + 4), c("vector", 
        "length", "mode", "content")))
    sumry[1, ] <- c("$lw", length(x$lw), mode(x$lw), "row weigths")
    sumry[2, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[3, ] <- c("$blo", length(x$blo), mode(x$blo), "column numbers")
    sumry[4, ] <- c("$tabw", length(x$tabw), mode(x$tabw), "array weights")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(3, 4), list((ntab + 5):(ntab + 7), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "Factors Table number Line number")
    sumry[2, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "Factors Table number Col number")
    sumry[3, ] <- c("$T4", nrow(x$T4), ncol(x$T4), "Factors Table number 1234")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    cat((ntab + 8), "$call: ")
    print(x$call)
    cat("\n")
    cat("names :\n")
    for (i in 1:ntab) {
        cat(names(x)[i], ":", names(x[[i]]), "\n")
    }
    cat("\n")
    indica <- as.factor(rep(1:ntab, x$blo))
    w <- split(x$cw, indica)
    cat("Col weigths :\n")
    for (i in 1:ntab) {
        cat(names(x)[i], ":", w[[i]], "\n")
    }
    cat("\n")
    cat("Row weigths :\n")
    cat(x$lw)
    cat("\n")
}
########### c.ktab" ########### 
"c.ktab" <- function (...) {
    x <- list(...)
    n <- length(x)
    if (any(lapply(x, class) != "ktab")) 
        stop("arguments imply object without 'ktab' class")
    nr <- unlist(lapply(x, function(x) nrow(x[[1]])))
    if (length(unique(nr)) != 1) 
        stop("arguments imply object with non constant row numbers")
    lw <- x[[1]]$lw
    nr <- length(lw)
    noms <- row.names(x[[1]][[1]])
    res <- NULL
    cw <- NULL
    blocks <- NULL
    for (i in 1:n) {
        if (any(x[[i]]$lw != lw)) 
            stop("arguments imply object with non constant row weights")
        if (any(row.names(x[[i]][[1]]) != noms)) 
            stop("arguments imply object with non constant row.names")
        blo.i <- x[[i]]$blo
        nblo.i <- length(blo.i)
        res <- c(res, unclass(x[[i]])[1:nblo.i])
        cw <- c(cw, x[[i]]$cw)
        blocks <- c(blocks, blo.i)
    }
    names(res) <- make.names(names(res), TRUE)
    res$lw <- lw
    res$cw <- cw
    res$blo <- blocks
    ktab.util.addfactor(res) <- list(blocks, length(lw))
    res$call <- match.call()
    class(res) <- "ktab"
    return(res)
}

########### t.ktab" ########### 
"t.ktab" <- function (x) {
    if (!inherits(x, "ktab")) 
        stop("object 'ktab' expected")
    blocks <- x$blo
    nblo <- length(blocks)
    res <- x
    r.n <- row.names(x[[1]])
    for (i in 1:nblo) {
        r.new <- row.names(x[[i]])
        if (any(r.new != r.n)) 
            stop("non equal row.names among array")
    }
    if (length(unique(blocks)) != 1) 
        stop("non equal col numbers among array")
    c.n <- names(x[[1]])
    for (i in 1:nblo) {
        c.new <- names(x[[i]])
        if (any(c.new != c.n)) 
            stop("non equal col.names among array")
    }
    new.row.names <- names(x[[1]])
    indica <- as.factor(rep(1:nblo, blocks))
    w <- split(x$cw, indica)
    col.w <- w[[1]]
    for (i in 1:nblo) {
        col.w.new <- w[[i]]
        if (any(col.w != col.w.new)) 
            stop("non equal column weights among array")
    }
    for (j in 1:nblo) {
        w <- x[[j]]
        w <- data.frame(t(w))
        row.names(w) <- new.row.names
        res[[j]] <- w
        blocks[j] <- ncol(w)
    }
    res$lw <- col.w
    res$cw <- rep(x$lw, nblo)
    res$blo <- blocks
    ktab.util.addfactor(res) <- list(blocks, length(res$lw))
    res$call <- match.call()
    class(res) <- "ktab"
    return(res)
}

########### row.names.ktab ########### 
"row.names.ktab" <- function (x) {
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    cha <- attr(x[[1]], "row.names")
    for (i in 1:ntab) {
        if (any(attr(x[[i]], "row.names") != cha)) 
            warnings(paste("array", i, "and array 1 have different row.names"))
    }
    return(cha)
}
########### row.names<-.ktab ########### 
"row.names<-.ktab" <- function (x, value) {
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    old <- attr(x[[1]], "row.names")
    if (!is.null(old) && length(value) != length(old)) 
        stop("invalid row.names length")
    value <- as.character(value)
    if (any(duplicated(value))) 
        stop("duplicate row.names are not allowed")
    for (i in 1:ntab) {
        attr(x[[i]], "row.names") <- value
    }
    x
}
########### col.names ########### 
"col.names" <- function (x) UseMethod("col.names")

########### col.names<- ########### 
"col.names<-" <- function (x, value) UseMethod("col.names<-")

########### col.names.ktab ########### 
"col.names.ktab" <- function (x) {
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    cha <- unlist(lapply(1:ntab, function(y) attr(x[[y]], "names")))
    return(cha)
}
########### col.names<-.ktab ########### 
"col.names<-.ktab" <- function (x, value) {
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    old <- unlist(lapply(1:ntab, function(y) attr(x[[y]], "names")))
    if (!is.null(old) && length(value) != length(old)) 
        stop("invalid col.names length")
    value <- as.character(value)
    indica <- as.factor(rep(1:ntab, x$blo))
    for (i in 1:ntab) {
        if (any(duplicated(value[indica == i]))) 
            stop("duplicate col.names are not allowed in the same array")
        attr(x[[i]], "names") <- value[indica == i]
    }
    x
}


########### tab.names ########### 
# fonction g�n�rique
"tab.names" <- function (x) UseMethod("tab.names")
########### tab.names.ktab ########### 
# m�thode pour ktab
"tab.names.ktab" <- function (x) {
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    cha <- names(x)[1:ntab]
    return(cha)
}
########### tab.names<- ########### 
# fonction g�n�rique
"tab.names<-" <- function (x, value) UseMethod("tab.names<-")
########### tab.names<-.ktab ########### 
# m�thode pour ktab
# les tab.names d'un ktab est le vecteur des noms des k premi�res composantes
# ce nombre de tableaux est la longueur de la composante blo
"tab.names<-.ktab" <- function (x, value) {
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    old <- tab.names(x)[1:ntab]
    if (!is.null(old) && length(value) != length(old)) 
        stop("invalid tab.names length")
    value <- as.character(value)
    if (any(duplicated(value))) 
        stop("duplicate tab.names are not allowed")
    names(x)[1:ntab] <- value
    x
}
########### ktab.util.names ###########
# utilitaire qui r�cup�re dans un ktab
# une liste de 3 �l�ments
# les noms des lignes "." les noms des tableaux
# les noms des colonnes sans duplicats
# les noms des tableaux "." 1234
# pour donner des �tiquettes aux TL, TC et T4 dans les graphiques
"ktab.util.names" <- function (x) {
    w <- row.names(x)
    w1 <- paste(w, as.character(x$TL[, 1]), sep = ".")
    w <- col.names(x)
    if (any(duplicated(w))) 
        w <- paste(w, as.character(x$TC[, 1]), sep = ".")
    w2 <- w
    w <- tab.names(x)
    l0 <- length(w)
    w3 <- paste(rep(w, rep(4, l0)), as.character(1:4), sep = ".")
    return(list(row = w1, col = w2, tab = w3))
}

########### ktab.util.addfactor<- ########### 
# utilitaire utilis� dans les ktab
# ajoute les composantes TL TC et T4
# x est un ktab presque achev�
# value est une liste contenant le vecteur des blocs de colonnes
# et le nombre de lignes
# on r�cup�re avec le nombre de tableaux, le nombre de variables par tableaux
# et le nombre de lignes en commun
"ktab.util.addfactor<-" <- function (x, value) {
    blocks <- value[[1]]
    nlig <- value[[2]]
    nblo <- length(blocks)
    w <- cbind.data.frame(gl(nblo, nlig), factor(rep(1:nlig, 
        nblo)))
    names(w) <- c("T", "L")
    x$TL <- w
    w <- NULL
    for (i in 1:nblo) w <- c(w, 1:blocks[i])
    w <- cbind.data.frame(factor(rep(1:nblo, blocks)), factor(w))
    names(w) <- c("T", "C")
    x$TC <- w
    w <- cbind.data.frame(gl(nblo, 4), factor(rep(1:4, nblo)))
    names(w) <- c("T", "4")
    x$T4 <- w
    x
}
