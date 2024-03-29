\name{corGrafen}
\alias{corGrafen}
\alias{coef.corGrafen}
\alias{corMatrix.corGrafen}
\title{Grafen's (1989) Correlation Structure}
\usage{
corGrafen(value, phy, form=~1, fixed = FALSE)
\method{coef}{corGrafen}(object, unconstrained = TRUE, ...)
\method{corMatrix}{corGrafen}(object,
                  covariate = getCovariate(object), corr = TRUE, ...)
}
\arguments{
  \item{value}{The \eqn{\alpha}{rho} parameter}
  \item{phy}{An object of class \code{phylo} representing the phylogeny
    (branch lengths are ignored) to consider}
  \item{object}{An (initialized) object of class \code{corGrafen}}
  \item{corr}{a logical value. If 'TRUE' the function returns the
    correlation matrix, otherwise it returns the variance/covariance
    matrix.}
  \item{fixed}{an optional logical value indicating whether the
    coefficients should be allowed to vary in the optimization, or kept
    fixed at their initial value. Defaults to 'FALSE', in which case the
    coefficients are allowed to vary.}
  \item{form}{ignored for now.}
  \item{covariate}{ignored for now.}
  \item{unconstrained}{a logical value. If 'TRUE' the coefficients are
    returned in unconstrained form (the same used in the optimization
    algorithm). If 'FALSE' the coefficients are returned in "natural",
    possibly constrained, form. Defaults to 'TRUE'}
  \item{...}{some methods for these generics require additional
    arguments. None are used in these methods.}
}
\description{
  Grafen's (1989) covariance structure. Branch lengths are computed using
  Grafen's method (see \code{\link{compute.brlen}}). The covariance
  matrice is then the traditional variance-covariance matrix for a
  phylogeny.
}
\value{
  An object of class \code{corGrafen} or the rho coefficient from an
  object of this class or the correlation matrix of an initialized
  object of this class.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
  \code{\link{corClasses}}, \code{\link{compute.brlen}}, \code{\link{vcv.phylo}}.
}
\references{
  Grafen, A. (1989) The phylogenetic regression. \emph{Philosophical
    Transactions of the Royal society of London. Series B. Biological
    Sciences}, \bold{326}, 119--157.
}
\keyword{models}
