\name{summary.pic.lm}
\alias{summary.pic.lm}

\title{ Summary function for pic.lm }
\description{
The usual \code{summary} does not work correctly with models produced by \code{pic.lm} because 
the hack by which the model degrees of freedom are reduced is detected and overruled.  This 
function is a slightly crippled version of \code{summary.lm}, altered to make it more trusting of 
objects passed to it. I think I probably need a new class of object really.}
\usage{
summary.pic.lm(object, correlation = FALSE, symbolic.cor = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"lm"}, that should be produced by \code{pic.lm}. }
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed}
  \item{symbolic.cor}{ logical. See \code{summary.lm} }
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
I'm too ashamed to give details.
}
\value{
As in \code{summary.lm}, but with d.f. reduced so there is only one per node in the phylogeny 
passed to \code{pic.lm}.
}

\author{Andy Purvis}

\section{Version }{1.1}
\section{Date }{16 January 2006}
\section{Modifications}{v1.1: Modified by Andy Purvis to correct failure to print coefficients 
table and to correct adjusted r-squared values}

\section{Warning }{If \code{summary.lm} is used, the wrong d.f. will be returned}


\seealso{ \code{\link[base]{summary.lm}}, \code{\link{pic.lm}} }


## The function is currently defined as
function (object, correlation = FALSE, symbolic.cor = FALSE,
    ...) 

#Slightly modified from R's base summary.lm function
#Modification sets df to that specified

{
    z <- object
    Qr <- object$qr
    if (is.null(z$terms) || is.null(Qr)) 
        stop("invalid 'lm' object:  no terms or qr component")
    n <- NROW(Qr$qr)
    p <- z$rank
    rdf <- z$df.residual
    if (rdf != z$df.residual) 
        warning("Degrees of freedom were modified")
    p1 <- 1:p
    r <- z$resid
    f <- z$fitted
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept")) 
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
        rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept")) 
            1
        else 0
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * (((rdf+length(tval)) - 
            df.int)/rdf) #Modified from summary.lm, which had n instead of rdf+length(tval)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
            numdf = p - df.int, dendf = rdf)
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
        1)]
    if (correlation) {
        ans$correlation <- (R * resvar)/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.lm"
    ans
  }
}
\keyword{ regression}
\keyword{ robust}
\keyword{ htest}