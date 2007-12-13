\name{subset-methods}
\docType{methods}
\alias{subset}
\alias{subset-methods}
\alias{subset,phylo-method}
\alias{subset,phylo4-method}
\alias{subset,phylo4d-method}
\title{Methods for creating subsets of phylogenies}
\description{
  Methods for creating subsets of phylogenies, based on pruning a tree to a set of terminal taxa, eventually to include pruning to include MRCA of taxa sets or based on a subtree rooted at a given node.
}
\section{Methods}{
  \describe{
    \item{phy = "phylo4d"}{drop tips and corresponding node and tip data}
    \item{phy = "phylo4"}{drop tips}
    \item{phy = "phylo"}{drop tips}
  }
}
\note{At the moment, this simply wraps \code{ape::drop.tip}.
    Renamed from \code{drop.tip} to \code{prune} in order to
    avoid conflicts with ape.}
\usage{\S4method{prune}{phylo4}(phy, tip, trim.internal = TRUE, subtree = FALSE,
         root.edge = 0)
       DropTip(phy,tip,\dots)
}
\arguments{
  \item{phy}{an object of class \code{"phylo4"}}
  \item{tip}{a vector of mode numeric or character specifying the tips
    to delete}
}
\details{
  To be described.

  The arguments \code{tip} can be either character or numeric. In the
  first case, it gives the labels of the tips to be deleted; in the
  second case the numbers of these labels in the vector
  \code{phy\$tip.label} are given

}
\value{
  an object of class \code{"phylo4"}
}
\author{Steven Kembel \email{skembel@berkeley.edu}}
\examples{

}
\keyword{manip}
\keyword{methods}

