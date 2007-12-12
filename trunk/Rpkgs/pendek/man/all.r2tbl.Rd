\name{all.r2tbl}
\alias{all.r2tbl}
\title{Compute distances between the root and each tip of a phylogeny}
\description{
Computes the total distance between the root of a phylogeny and each of its tips}
\usage{
all.r2tbl(phylogeny)
}

\arguments{
  \item{phylogeny}{A valid item of class \code{"phylo"}}
}
\details{Can be used instead of \code{\link{get.root.to.tip}}, and it is not
obvious to me why I wrote both.}
}
\value{
A vector in the same order as \code{phylogeny$tip.label}}

\author{Andy Purvis}

\section{Version }{1.3}
\section{Date }{19 May 2006}

\section{Warning }{Does not currently check that the phylogeny passed is valid}

\keyword{ manip }
