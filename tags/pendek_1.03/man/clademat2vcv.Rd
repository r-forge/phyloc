\name{clademat2vcv}
\alias{clademat2vcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Converts a clade matrix to VCV matrix.}
\description{
Takes a 'clade.matrix' object and converts it into a variance covariance matrix.}
\usage{
clademat2vcv(clmat, method = "loop")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clmat}{ An object of class 'clade.matrix'.}
  \item{method}{ One of either 'loop' or 'bycolumn', with 'bycolumn' being faster. }
}
\details{
An experimental function, written in an attempt to get a quicker way than \code{\link{vcv.phylo}} of turning large phylogenies into VCV matrices. The 'method' arguments reflects two possible solutions. The user should beware that this  has not been fully tested. 

}
\value{
A matrix showing the the variance (diagonal) and covariance (off-diagonal) between tips in a phylogeny with branch lengths.
}
\author{David Orme}

\section{Warning}{Not fully tested}

\seealso{ \code{\link{clade.matrix}}, \code{\link{vcv.phylo}} }
\examples{
data(bird.orders)
bird.cm <- clade.matrix(bird.orders)
clademat2vcv(bird.cm)
}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{utilities}% __ONLY ONE__ keyword per line
