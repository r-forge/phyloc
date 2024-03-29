\name{theta.s}
\alias{theta.s}
\title{Population Parameter THETA using Segregating Sites
  in DNA Sequences}
\usage{
theta.s(s, n, variance = FALSE)
}
\arguments{
  \item{s}{a numeric giving the number of segregating sites.}
  \item{n}{a numeric giving the number of sequences.}
  \item{variance}{a logical indicating whether the variance of the
    estimated THETA should be returned (\code{TRUE}), the default being
    \code{FALSE}.}
}
\description{
  This function computes the population parameter THETA using the
  number of segregating sites \code{s} in a sample of \code{n} DNA sequences.
}
\value{
  a numeric vector of length one with the estimated theta (the default),
  or of length two if the standard error is returned
  (\code{variance = TRUE}).
}
\note{
  The number of segregating sites needs to be computed beforehand, for
  instance with the function \code{seg.sites} (see example below).
}
\references{
  Watterson, G. (1975) On the number of segragating sites in genetical
  models without recombination. \emph{Theoretical Population Biology},
  \bold{7}, 256--276.

  Tajima, F. (1989) Statistical method for testing the neutral mutation
  hypothesis by DNA polymorphism. \emph{Genetics}, \bold{123}, 585--595.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{theta.h}}, \code{\link{theta.k}}, \code{\link{seg.sites}},
  \code{\link{nuc.div}}
}
\examples{
data(woodmouse)
y <- seg.sites(woodmouse)
s <- length(y)
n <- length(woodmouse)
theta.s(s, n)
theta.s(s, n, variance = TRUE)
}
\keyword{manip}
\keyword{univar}
