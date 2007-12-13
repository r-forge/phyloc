\name{coerce-methods}
\docType{methods}
\alias{coerce-methods}
\alias{coerce,phylo,phylo4-method}
\alias{coerce,multiPhylo4,multiPhylo-method}
\alias{coerce,multiPhylo,multiPhylo4-method}
\alias{coerce,phylo4,phylo-method}
\title{Converting between phylo4 and ape trees}
\description{
  Conversion functions to convert betwee \code{ape} objects
  (\code{phylo}, \code{multiPhylo}) and \code{phylo4} objects
  (\code{phylo4}, \code{phylo4d}, \code{multiPhylo4d})}
\section{Methods}{
  \describe{
    \item{from = "phylo", to = "phylo4"}{ ~~describe this method here }
    \item{from = "multiPhylo4", to = "multiPhylo"}{ ~~describe this method here }
    \item{from = "multiPhylo", to = "multiPhylo4"}{ ~~describe this method here }
    \item{from = "phylo4", to = "phylo"}{ ~~describe this method here }
  }
}
\author{Ben Bolker}
\examples{
    library(ape)
  example(read.tree)
## round trip conversion
  p1 = tree.owls
  P1 = as(tree.owls,"phylo4")
  P1
  identical(p1,as(P1,"phylo"))
}
\keyword{methods}
