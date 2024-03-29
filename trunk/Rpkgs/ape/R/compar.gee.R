## compar.gee.R (2006-10-11)

##   Comparative Analysis with GEEs

## Copyright 2002-2006 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

compar.gee <- function(formula, data = NULL, family = "gaussian", phy,
                       scale.fix = FALSE, scale.value = 1)
{
    if (is.null(data)){ 
        data <- parent.frame() 
        for(i in attr(terms.formula(formula), 'term.labels')) {
            if(!.ape.quiet && !hasNames(get(i))) {
                warning(paste("Formula argument '", i, "' does not have names. Order is assumed to the be the same as phy$tip.label", sep = ""))
            }
        }
    }
    else {
        if(!any(is.na(match(rownames(data), phy$tip.label))))
          data <- data[phy$tip.label, ]
        else warning("the rownames of the data.frame and the names of the tip labels
do not match: the former were ignored in the analysis.")
    }
    effect.assign <- attr(model.matrix(formula, data = data), "assign")
    for (i in all.vars(formula)) {
        if (any(is.na(eval(parse(text = i), envir = data))))
          stop("the present method cannot (yet) be used directly with missing data: you may consider removing the species with missing data from your tree with the function `drop.tip'.")
    }
    if (is.null(phy$edge.length))
      stop("the tree has no branch lengths.")
    R <- vcv.phylo(phy, cor = TRUE)
    id <- rep(1, dim(R)[1])
    geemod <- do.call("gee", list(formula, id, data = data, family = family, R = R,
                                  corstr = "fixed", scale.fix = scale.fix,
                                  scale.value = scale.value))
    W <- geemod$naive.variance
    if (family == "binomial")
      W <- summary(glm(formula, family = quasibinomial, data = data))$cov.scaled
    N <- geemod$nobs
    dfP <- sum(phy$edge.length)*N / sum(diag(vcv.phylo(phy)))
    obj <- list(call = geemod$call,
                effect.assign = effect.assign,
                nobs = N,
                coefficients = geemod$coefficients,
                residuals = geemod$residuals,
                family = geemod$family$family,
                link = geemod$family$link,
                scale = geemod$scale,
                W = W,
                dfP = dfP)
    class(obj) <- "compar.gee"
    obj
}

print.compar.gee <- function(x, ...)
{
    nas <- is.na(x$coef)
    coef <- x$coef[!nas]
    cnames <- names(coef)
    coef <- matrix(rep(coef, 4), ncol = 4)
    dimnames(coef) <- list(cnames,
                           c("Estimate", "S.E.", "t", "Pr(T > |t|)"))
    df <- x$dfP - dim(coef)[1]
    coef[, 2] <- sqrt(diag(x$W))
    coef[, 3] <- coef[, 1]/coef[, 2]
    if (df < 0) {
        warning("not enough degrees of freedom to compute P-values.")
        coef[, 4] <- NA
    } else coef[, 4] <- 2 * (1 -  pt(abs(coef[, 3]), df))
    residu <- quantile(as.vector(x$residuals))
    names(residu) <- c("Min", "1Q", "Median", "3Q", "Max")
    cat("\nCall:\n")
    cat("  formula: ")
    print(x$call$formula)
    cat("\nNumber of observations: ", x$nobs, "\n")
    cat("\nModel:\n")
    cat(" Link:                     ", x$link, "\n")
    cat(" Variance to Mean Relation:", x$family, "\n")
    cat("\nSummary of Residuals:\n")
    print(residu)
    if (any(nas))
        cat("\n\nCoefficients: (", sum(nas), " not defined because of singularities)\n",
            sep = "")
    else cat("\n\nCoefficients:\n")
    print(coef)
    cat("\nEstimated Scale Parameter: ", x$scale)
    cat("\n\"Phylogenetic\" df (dfP): ", x$dfP, "\n")
}

drop1.compar.gee <- function(object, scope, quiet = FALSE, ...)
{
    fm <- formula(object$call)
    trm <- terms(fm)
    z <- attr(trm, "term.labels")
    ind <- object$effect.assign
    n <- length(z)
    ans <- matrix(NA, n, 3)
    for (i in 1:n) {
        wh <- which(ind == i)
        ans[i, 1] <- length(wh)
        ans[i, 2] <- t(object$coefficients[wh]) %*%
          solve(object$W[wh, wh]) %*% object$coefficients[wh]
    }
    df <- object$dfP - length(object$coefficients)
    if (df < 0) warning("not enough degrees of freedom to compute P-values.")
    else ans[, 3] <- pf(ans[, 2], ans[, 1], df, lower.tail = FALSE)
    colnames(ans) <- c("df", "F", "Pr(>F)")
    rownames(ans) <- z
    if (any(attr(trm, "order") > 1) && !quiet)
      warning("there is at least one interaction term in your model:
you should be careful when interpreting the significance of the main effects.")
    class(ans) <- "anova"
    attr(ans, "heading") <- c("Single term deletions\n\nModel:\n",
                              as.character(as.expression(fm)))
    ans
}
