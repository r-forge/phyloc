\name{mrca}
\alias{mrca}
\title{Find Most Recent Common Ancestors Between Pairs}
\description{
  This function returns for each pair of tips (and nodes) its most
  recent common ancestor (MRCA).
}
\usage{
mrca(phy, full = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{full}{a logical indicating whether to return the MRCAs among
    all tips and nodes (if \code{TRUE}); the default is to return only
    the MRCAs among tips.}
}
\details{
  The diagonal is set to the number of the tips (and nodes if \code{full
    = TRUE}).

  If \code{full = FALSE}, the colnames and rownames are set with the tip
  labels of the tree; otherwise the numbers are given as names.
}
\value{
  a matrix of mode numeric.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\keyword{manip}
