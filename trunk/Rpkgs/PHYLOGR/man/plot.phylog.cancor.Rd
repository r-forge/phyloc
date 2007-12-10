\name{plot.phylog.cancor}
\alias{plot.phylog.cancor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot a phylog.cancor object }
\description{
  Plots histogram of the canonical correlations for simulated data
  as returned from a phylog.cancor object; with
  vertical bars indicates the values from the original (''real'') data
  (the one with sim.counter=0), and in parenthesis their
  'correlation-wise' p-value (see summary.phylog.cancor).

}
\usage{
\method{plot}{phylog.cancor}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ an object of class phylog.cancor returned from a previous call to
 cancor.phylog. }
 \item{...}{ other parameters to be passed to through to plotting
   functions (currently not used). }

}

\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.

}
\author{ Ramon Diaz-Uriarte and Theodore Garland, Jr.}
\section{WARNING}{These histograms are in the spirit of the 'correlation-wise'
  p-values returned from summary.phylog.cancor; see Details for summary.phylog.cancor.}
\seealso{\code{\link{cancor.phylog}}, \code{\link{summary.phylog.cancor}}}

\examples{
data(SimulExample)
ex1.cancor <- cancor.phylog(SimulExample[,c(1,2,3,4,5)],SimulExample[,c(1,2,6,7,8)])
ex1.cancor
summary(ex1.cancor)
par(mfrow=c(1,3))
plot(ex1.cancor)
}
\keyword{methods}%-- one or more ...

