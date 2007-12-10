\name{plot.phylog.lm}
\alias{plot.phylog.lm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot a phylog.lm object }
\description{
  Plots histogram of the the fitted coefficients and F-tests for simulated data
  as returned from a phylog.lm object; with
  vertical bars indicates the values from the original (''real'') data
  (the one with sim.counter=0). For F-values, the number in parenthesis
  indicates the p-value (see summary.phylog.lm).
  

}
\usage{
\method{plot}{phylog.lm}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ an object of class phylog.lm returned from a previous call to
   lm.phylog. }
 \item{...}{ other parameters to be passed to through to plotting
   functions (currently not used). }

}

\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
}
\author{ Ramon Diaz-Uriarte and Theodore Garland, Jr.}

\seealso{\code{\link{lm.phylog}}, \code{\link{summary.phylog.lm}}}

\examples{
data(SimulExample)
ex1.lm <- lm.phylog(y ~ x1 + diet, weights=x2, max.num=20, data=SimulExample)
options(digits=5)
plot(ex1.lm)
}
\keyword{methods}%-- one or more ...

