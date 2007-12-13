\name{simulate.brownian}
\alias{simulate.brownian}
\title{Simulation of character evolution by Brownian motion}
\description{
This function simulates the evolution of a single quantitative character
by Brownian motion along the branches of the phylogeny passed to it.
}
\usage{
simulate.brownian(phylogeny)
}
\arguments{
  \item{phylogeny}{An item of class \code{"phylo"}}
}
\details{
The amount and sign of change in the character is computed for each edge in the phylogeny.
A vector of the changes is added to the phylogeny object, which is then returned.  This function is a special case of the more general
\code{\link{simulate.with.kappa}}, with \code{kappa} = 1.
}
\value{
An object of class \code{"phylo"} with an extra vector, \code{"changes"}.
}
\references{Felsenstein, J. 1985  \emph{"American Naturalist"} \bold{"125"}, 1--15.}
\author{Andy Purvis}

\section{Warning }{At present, this function does not check that it has been passed a valid 
object of class \code{"phylo"}.}

\seealso{\code{\link{simulate.with.kappa}}}

## The function is currently defined as
function(phylogeny){
#Add a vector of changes simulated under a Brownian motion model
rands<-rnorm(length(phylogeny$edge.length))
phylogeny$changes<-rands*sqrt(phylogeny$edge.length)
return(phylogeny)
  }
}
\keyword{manip}