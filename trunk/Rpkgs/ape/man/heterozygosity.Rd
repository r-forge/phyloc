\name{heterozygosity}
\alias{heterozygosity}
\alias{H}
\title{Heterozygosity at a Locus Using Gene Frequencies}
\usage{
heterozygosity(x, variance = FALSE)
H(x, variance = FALSE)
}
\arguments{
  \item{x}{a vector or a factor.}
  \item{variance}{a logical indicating whether the variance of the 
    estimated heterozygosity should be returned (\code{TRUE}), the
    default being \code{FALSE}.}
}
\description{
  This function computes the mean heterozygosity from gene frequencies,
  and returns optionally the associated variance.
}
\value{
  a numeric vector of length one with the estimated mean heterozygosity
  (the default), or of length two if the variance is returned
  \code{variance = TRUE}.
}
\details{
  The argument \code{x} can be either a factor or a vector. If it is a
  factor, then it is taken to give the individual alleles in the
  population. If it is a numeric vector, then its values are taken to be
  the numbers of each allele in the population. If it is a non-numeric
  vector, it is a coerced as a factor.

  The mean heterozygosity is estimated with:

  \deqn{\hat{H} = \frac{n}{n-1} \left(1 - \sum_{i=1}^k p_i^2 \right)}{%
    H = n(1 - SUM (FROM i=1 TO k) p_i^2)/(n - 1)}

  where \eqn{n} is the number of genes in the sample, \eqn{k} is the
  number of alleles, and \eqn{p_i} is the observed (relative) frequency
  of the allele \eqn{i}.
}
\references{
  Nei, M. (1987) \emph{Molecular evolutionary genetics}. New York:
  Columbia University Press.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{theta.s}}
}
\keyword{manip}
\keyword{univar}
