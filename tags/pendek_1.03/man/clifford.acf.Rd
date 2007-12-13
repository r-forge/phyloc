\name{clifford.acf}
\alias{clifford.acf}
\title{Calculate Autocorrelation Function}
\description{
Calculates a spatial autocorrelation function for a matrix using all lags within the matrix, based on the spatially corrected correlation test described by Clifford et al (1989).}
\usage{
clifford.acf(x, ew.wrap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric matrix.}
  \item{ew.wrap}{A logical indicating whether to allow the autocorrelation function to wrap lags across the left and right edge of the matrix.}
}
\details{
The function computes an acf matrix as the average product of matrix elements at each lag within the matrix. The \code{ew.wrap} wraps lags around the left and right edge of the matrix where elements are closer together across that edge than within the matrix. This is intended to more accurately account for spatial autocorrelation at a global level.

The number of rows and columns in the acf is (2n)-1, representing the lags from each focal element in the input matrix to each other element.}
\value{
 A list containing two components:
  \item{acf}{The autocorrelation function}
  \item{nok}{The number of lags contributing to each element in the \code{acf}.}
}

\references{Clifford P, Richardson S and D Hemon (1989) 'Assessing the siginificance of the correlation between two spatial processes.' Biometrics 45(1) 123 - 134.}
\author{David Orme}
\note{Ported from a Fortran implementation written by Jack Lennon.}

\seealso{ \code{\link{clifford.test}} }
\examples{
A <- matrix(rep(1:4,each=3), ncol=4) + rnorm(12)
clifford.acf(A)
clifford.acf(A, ew.wrap=TRUE)
}
\keyword{ manip }% at least one, from doc/KEYWORDS
\keyword{ spatial }% __ONLY ONE__ keyword per line
