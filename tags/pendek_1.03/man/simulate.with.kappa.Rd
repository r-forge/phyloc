\name{simulate.with.kappa}
\alias{simulate.with.kappa}
\title{Simulate character evolution using modified Brownian motion model}
\description{
Simulates character evolution in which variance is proportional to branch length
raised to a power.}
\usage{
simulate.with.kappa(phylogeny, kappa = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{a valid object of class \code{"phylo"}.}
  \item{kappa}{numeric. A value of 0 makes branch lengths equal; 1 (the default) is equivalent
  to Brownian motion.}
}
}
\value{
An object of class \code{"phylo"}, that includes an additional vector \code{changes} containing the
simulated changes in the character along each edge of the phylogeny.}
\references{Garland, T. Jr., Harvey, P. H., and Ives, A.R. 1992 \emph{Systematic Biology} \bold{41}, 18--32; Pagel, M. 1997 \emph{Zoologica Scripta} \bold{26}, 331-348}
\author{Andy Purvis}

\section{Warning }{Does not presently check that the arguments passed are valid.}

\seealso{\code{\link{simulate.brownian}}}

## The function is currently defined as
function(phylogeny,kappa=1){
#Add a vector of changes simulated under a modified Brownian model in which
#  branch lengths are raised to a power, kappa.
rands<-rnorm(length(phylogeny$edge.length))
phylogeny$changes<-rands*sqrt(phylogeny$edge.length^kappa)
return(phylogeny)
  }
}
\keyword{manip}