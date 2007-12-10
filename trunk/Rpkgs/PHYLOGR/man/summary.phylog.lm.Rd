\name{summary.phylog.lm}
\alias{summary.phylog.lm}

\title{ Summarize a phylog.lm object}
\description{
  A 'method' for objects of class phylog.lm. Summarizes the results of
  the fitted linear models for the simulated (and original) data and
  returns p-values. There is a print 'method' for this summary.
}
\usage{
\method{summary}{phylog.lm}(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{object}{ an object of class phylog.lm returned from a previous call to
   lm.phylog.}
 \item{...}{ further arguments passed to or from other methods
   (currently not used).}
}
\details{ The p-value is computed as the (number of simulated data sets
  with F-value larger than original (''real'') data + 1) / (number of simulated
  data sets + 1).
  The quantiles are obtained with the function quantile, and thus linear interpolation is used.

}
\value{ A list (of class summary.phylog.lm) with components
  \item{call}{the call to lm.phylog}
  \item{original.model}{the model fitted to the original data.}
  \item{original.Fvalue}{the F-values from the original data.}
  \item{p.value}{the p-values for the marginal F-tests and overall F.}
  \item{quant.F.value}{the quantiles of the distribution of F-values
    from the simulated data.}
  \item{num.simul}{the number of simulations used in the analyses}
}
\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
}
\author{ Ramon Diaz-Uriarte and Theodore Garland, Jr.}

\seealso{\code{\link{lm.phylog}}, \code{\link{plot.phylog.lm}}}

\examples{
data(SimulExample)
ex1.lm <- lm.phylog(y ~ x1 + diet, weights=x2, max.num=20, data=SimulExample)
summary(ex1.lm)
}
\keyword{methods}%-- one or more ...
