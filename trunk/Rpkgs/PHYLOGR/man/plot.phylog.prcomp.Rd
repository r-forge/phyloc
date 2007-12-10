\name{plot.phylog.prcomp}
\alias{plot.phylog.prcomp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot a phylog.prcomp object }
\description{
  Plots histogram of the eigenvalues for simulated data
  as returned from a phylog.prcomp object; with
  vertical bars indicates the values from the original (''real'') data
  (the one with sim.counter=0).

}
\usage{
\method{plot}{phylog.prcomp}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ an object of class phylog.prcomp such as can be obtained
   from a previous call to prcomp.phylog. }
 \item{...}{ other parameters to be passed to through to plotting
   functions (currently not used). }
}

\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.

}
\author{ Ramon Diaz-Uriarte and Theodore Garland, Jr.}
\section{WARNING}{These histograms are in the spirit of the 'naive'
  p-values returned from summary.phylog.prcomp; see Details for summary.phylog.prcomp}

\seealso{\code{\link{prcomp.phylog}}, \code{\link{summary.phylog.prcomp}}}

\examples{

data(SimulExample)
ex1.prcomp <- prcomp.phylog(SimulExample[,-11]) # 11th column is a factor
options(digits=5)
plot(ex1.prcomp)
}
\keyword{methods}%-- one or more ...

