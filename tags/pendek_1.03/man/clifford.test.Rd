\name{clifford.test}
\alias{clifford.test}
\alias{print.clifford.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clifford test for spatial correlation.}
\description{
Calculates a estimate of the effective sample size accounting for the degree of spatial autocorrelation between two matrices. This effective sample size is then used as the basis for a t-test of the significance of the Pearson correlation coefficient between the two variables.
}
\usage{
clifford.test(A, B, ew.wrap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{A numeric matrix.}
  \item{B}{A second numeric matrix with the same number of rows and columns as \code{A}. Ideally any NA values should match to NA values in \code{A}; if this isn't the case then the pairwise complete cases are used.}
  \item{ew.wrap}{A logical indicating whether to allow the autocorrelation function to wrap lags across the left and right edge of the matrix.}
}
\details{
The function implements the test described by Clifford et al. (1989) which describes a test of association based on the correlation coefficient between two spatially autocorrelated processes using data in a lattice. The test considers all lags within the data matrices to alter the effective sample size based on an estimate of standardized covariance. With positive autocorrelation, this is typically a reduction in sample size, but negative autocorrelation can give rise to an effective sample size that is larger than the number of non NA cases. The implementation of the autocorrelation function is not particularly quick!

The \code{ew.wrap} option is an extension that allows lags to wrap around from the left to the right hand edge of the matrix. This is intended to allow for the use of the test on global datasets.
}
\value{
A list of class 'clifford.test' containing the components:
  \item{A}{A list giving the name of matrix A, the matrix itself, the calculated autocorrelation function and the variance estimate.}
  \item{B}{A list giving for matrix B as above.}
  \item{nok}{A matrix giving the number of valid contributions at each lag.}
  \item{rpearson}{The Pearson correlation coefficient between the two matrices.}
  \item{n}{The number of non NA cases in the matrices.}
  \item{ncc}{The sum of the products of the \code{nok} and the two autocorrelation functions.}
  \item{var.xy}{The conditional variance between the two matrices.}
  \item{var.r}{The variance of the Pearson correlation coefficient.}
  \item{ess}{The effective sample size.}
  \item{w}{The standardized covariance.}
  \item{t}{The t statistic of \code{rpearson} given the effective sample size.}
  \item{p}{The significance of the correlation between the two matrices, accounting for the spatial autocorrelation.}
}
\references{Clifford P, Richardson S and D Hemon (1989) 'Assessing the siginificance of the correlation between two spatial processes.' Biometrics 45(1) 123 - 134.}
\author{David Orme}
\note{Ported from a Fortran implementation written by Jack Lennon.}


\seealso{\code{\link{clifford.acf}}}
\examples{
## Two matrices with spatial autocorrelation
A <- matrix(rep(1:10,each=5), ncol=10) + rnorm(50)
B <- matrix(rep(1:10,each=5), ncol=10) + rnorm(50)
clifford.test(A, B)
clifford.test(A, B, ew.wrap=TRUE)
}
\keyword{htest}% at least one, from doc/KEYWORDS
\keyword{spatial}% at least one, from doc/KEYWORDS
