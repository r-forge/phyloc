\name{seg.sites}
\alias{seg.sites}
\title{
  Find Segregating Sites in DNA Sequences
}
\usage{
seg.sites(x)
}
\arguments{
  \item{x}{a matrix or a list which contains the DNA sequences.}
}
\description{
  This function gives the indices of segregating (polymorphic) sites in
  a sample of DNA sequences.
}
\details{
  If the sequences are in a list, all the sequences must be of the same
  length.
}
\value{
  A numeric vector giving the indices of the segregating sites.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\note{
  The present version looks for the sites which are ``variable'' in the
  data in terms of different \emph{letters}. This may give unexpected
  results if there are ambiguous bases in the data.
}
\seealso{
  \code{\link{base.freq}}, \code{\link{GC.content}},
  \code{\link{theta.s}}, \code{\link{nuc.div}}
}
\examples{
data(woodmouse)
y <- seg.sites(woodmouse)
y
length(y)
}
\keyword{univar}
