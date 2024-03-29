\name{corMartins}
\alias{corMartins}
\alias{coef.corMartins}
\alias{corMatrix.corMartins}
\title{Martins's (1997) Correlation Structure}
\usage{
corMartins(value, phy, form=~1, fixed = FALSE)
\method{coef}{corMartins}(object, unconstrained = TRUE, ...)
\method{corMatrix}{corMartins}(object,
		covariate = getCovariate(object), corr = TRUE, ...)
}
\arguments{
  \item{value}{The \eqn{\alpha}{alpha} parameter}
  \item{phy}{An object of class \code{phylo} representing the phylogeny
    (with branch lengths) to consider}
  \item{object}{An (initialized) object of class \code{corMartins}}
  \item{corr}{a logical value. If 'TRUE' the function returns the
    correlation matrix, otherwise it returns  the variance/covariance
    matrix.}
  \item{fixed}{an optional logical value indicating whether the
    coefficients should be allowed to vary in the optimization, ok kept
    fixed at their initial value. Defaults to 'FALSE', in which case the
    coefficients are allowed to vary.}
  \item{form}{ignored for now.}
  \item{covariate}{ignored for now.}
  \item{unconstrained}{a logical value. If 'TRUE' the coefficients are returned
    in unconstrained form (the same used in the optimization
    algorithm). If 'FALSE' the coefficients are returned in
    "natural", possibly constrained, form. Defaults to 'TRUE'}
	\item{...}{some methods for these generics require additional arguments.
		None are used in these methods.}
}
\description{
	Martins and Hansen's (1997) covariance structure:
		\deqn{V_{ij} = \gamma \times e^{-\alpha t_{ij}}}{%
					Vij = gamma . exp(-alpha . tij)}
	where \eqn{t_{ij}}{tij} is the phylogenetic distance between taxa \eqn{i}{i} and \eqn{j}{j} and \eqn{\gamma}{gamma} is a constant.}
}
\value{
	An object of class \code{corMartins} or the alpha coefficient from an object of this class
	or the correlation matrix of an initialized object of this class.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
  \code{\link{corClasses}}.
}
\references{
  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the comparative
  method: a general approach to incorporating phylogenetic information
  into the analysis of interspecific data. \emph{American Naturalist},
  \bold{149}, 646--667.
}
\keyword{models}
