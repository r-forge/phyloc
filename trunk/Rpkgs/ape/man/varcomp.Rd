\name{varcomp}
\alias{varcomp}
\title{Compute Variance Component Estimates}
\description{
  Get variance component estimates from a fitted \code{lme} object.
}
\usage{
varcomp(x, scale = FALSE, cum = FALSE)
}
\arguments{
  \item{x}{A fitted \code{lme} object}
  \item{scale}{Scale all variance so that they sum to 1}
  \item{cum}{Send cumulative variance components.}
}
\details{
  Variance computations is done as in Venables and Ripley's book.
}
\value{
  A named vector of class \code{varcomp} with estimated variance components.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002) \emph{Modern applied statistics
  with S (fourth edition)}. New York: Springer-Verlag.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{ \code{\link[nlme]{lme}}, \code{\link{Correlogram}}}
\examples{
data(carnivora)
m <- lme(log10(SW) ~ 1, random = ~ 1|Order/SuperFamily/Family/Genus, data=carnivora)
v <- varcomp(m, TRUE, TRUE)
plot(v)
}
\keyword{regression}
\keyword{dplot}
