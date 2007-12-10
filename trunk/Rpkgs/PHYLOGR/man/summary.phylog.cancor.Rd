\name{summary.phylog.cancor}
\alias{summary.phylog.cancor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Summarize a phylog.cancor object}
\description{
  A 'method' for objects of class phylog.cancor. Shows the original
  data, and provides p-values and quantiles of the canonical
  correlations based on the simulated data. There is a print 'method'
  for this summary.
}
\usage{
\method{summary}{phylog.cancor}(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{object}{ an object of class phylog.cancor returned from a previous
   call to cancor.phylog. }
 \item{...}{further arguments passed to or from other methods (currently
   not used)}.
}
\details{
  To test the hypothesis that all population canonical correlations are
  zero we use the likelihood-ratio statistic from Krzanowski (pp. 447
  and ff.); this
  statistic is computed for the original data set and for each of the
  simulated data sets, and we obtain the p-value as  (number of simulated data sets
  with LR statistic larger than original (''real'') data + 1) / (number of simulated
  data sets + 1).  Note that a test of this same hypothesis using the
  Union-Intersection approach is equivalent to the test we implement
  below for the first canonical correlation.
  
  The p-values for the individual canonical correlations are calculated in two different ways. For the 'component-wise'
  ones the p-value for a particular correlation is
  (number of simulated data sets
  with canonical correlation larger than original (''real'') data + 1) / (number of simulated
  data sets + 1). With this approach, you can find that the p-value for,
  say, the second canonical correlation is smaller than the first, which
  is not sensible. It only makes sense to examine the second
  canonical correlation if the first one is ''significant'', etc. Thus,
  when considering the significance of the second canonical correlation
  we should account for the value of the first. In other words, there is only
  support against the null hypothesis (of no singificant second
  canonical correlation) if both the first and the second canonical
  correlations from the observed data set are larger than most of the
  simulated data sets. We can account for what happens with
  the first canonical correlation by computing the p-value of the second
  canonical correlation as the number of simulations in which the second
  simulated canonical correlation is larger than the observed, or the
  first simulated canonical correlation is larger than the observed one,
  or both (so that the only cases that count agains the null are those
  where both the first ans second canonical correlations are smaller
  than the observed ones); these we call 'Multiple' p-values.

}

\value{ A list (of class summary.phylog,cancor) with elements
  \item{call}{the call to function cancor.phylog.}
  \item{original.LR.statistic} {the likelihood ratio statistic for the
    test that all canonical correlations are zero}
  \item{original.canonicalcorrelations}{the canonical correlations
    corresponding to the original (''real'') data set.}
  \item{p.value.overall.test} {the p-value for the test that all
    canonical correlations are zero}
  \item{p.value.corwise}{the correlation-wise p-value ---see Details}
  \item{p.value.mult}{the multiple correlations p-value; see Details}
  \item{quant.canonicalcorrelations}{the quantiles from the simulated
    canonical correlations; linear interpolation is used. Note that
    these quantiles are in the spirit of the  ''naive p.values''.}
  \item{num.simul}{the number of simulations used in the analyses}
}


\references{

  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
  
  Krzanowski, W. J. (1990) \emph{Principles of multivariate analysis}
  Oxford University Press.
  
  Manly,B. F. J. (1997) \emph{Randomization, bootstraping, and Monte Carlo
    methods in biology, 2nd ed.} Chapman & Hall.
  
  Morrison, D. F. (1990) \emph{Multivariate statistcal methods, 3rd ed.} McGraw-Hill.
}
\author{Ramon Diaz-Uriarte and Theodore Garland, Jr.}

\section{WARNING}{
  It is
  necessary to be careful with the null hypothesis you are testing and
  how the null data set ---the simulations--- are generated. For
  instance, suppose you want to examine the canonical correlations
  between sets x and y; you will probably want to generate x and y each
  with the observed correlations within each set so that the
  correlations within each set are maintained (but with no
  correlations among sets). You probably do not want to generate each of
  the x's as if they were independent of each other x, and ditto for y, since
  that will destroy the correlations within each set; see some
  discussion in Manly, 1997.

}

\seealso{ \code{\link{read.sim.data}},  \code{\link{summary.phylog.cancor}}}

\examples{
data(SimulExample)
ex1.cancor <- cancor.phylog(SimulExample[,c(1,2,3,4,5,6)],SimulExample[,c(1,2,7,8)])
summary(ex1.cancor)

}
\keyword{multivariate}%-- one or more ...
\keyword{methods}

