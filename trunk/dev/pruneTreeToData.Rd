\name{pruneTreeToData}
\alias{pruneTreeToData}
\title{Prune a phylo tree to include only species in a data object}
\usage{pruneTreeToData(x,tree,NA.omit=TRUE,...)}
\arguments{
  \item{x}{a vector, data.frame or matrix object}
  \item{tree}{a phylo object}
  \item{NA.omit}{should species present in the data but with NA values be included? (TRUE or FALSE)}
}
\description{
  \code{pruneTreeToData} is a utility function for pruning a phylo tree to include only the taxa present in the supplied data object. If the data are in a vector, the names of the vector are matched, while for a data.frame or matrix the row.names are matched.
}
\value{
 A tree pruned to include only terminal taxa present in the supplied data object.
}
\author{Steven Kembel \email{skembel@berkeley.edu}}

\keyword{manip}
