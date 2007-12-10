\name{cancor.phylog}
\alias{cancor.phylog}
\title{Canonical Correlation Analysis from Simulated Data Sets}
\description{
  Performs a canonical correlation analysis on two sets of simulated data
  and returns the canonical correlations.
}
\usage{
cancor.phylog(data1, data2, max.num=0, exclude.tips=NULL,lapply.size=100,
              xcenter=TRUE, ycenter=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data1}{the columns from a data set returned from read.sim.data that you want to use as the first set in the canonical correlation analysis. The first column MUST be sim.counter and the second Tips.}
 \item{data2}{ the columns from a data set returned from read.sim.data that you want to use as the second set in the canonical correlation analysis. The first column MUST be sim.counter and the second Tips.}
  \item{max.num}{ if different from 0, the maximum number of simulations to analyze. }
 \item{exclude.tips}{ an optional vector giving the names of tips to exclude from the analyses.}
 \item{lapply.size}{a tuning parameter that can affect the speed of calculations; see Details in phylog.lm.}
 \item{xcenter}{ should the x columns be centered? defaults to yes; see cancor.}
 \item{ycenter}{ should the y columns be centered? defaults to yes; see
   cancor.}
}


\value{ A list (of class phylog.cancor) with components
  \item{call}{the call to the function}
  \item{CanonicalCorrelations}{the canonical correlations.  The one with
    sim.counter=0 corresponds to the original (''real'') data.}
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

\section{WARNING}{It is
  necessary to be careful with the null hypothesis you are testing and
  how the null data set ---the simulations--- are generated. For
  instance, suppose you want to examine the canonical correlations
  between sets x and y; you will probably want to generate x and y each
  with the observed correlations within each set so that the
  correlations within each set are maintained (but with no
  correlations among sets). You probably do not want to generate each of
  the x's as if they were independent of each other x, and ditto for y, since
  that will destroy the correlations within each set; see some
  discussion in Manly, 1997.}

\seealso{ \code{\link{read.sim.data}},  \code{\link{summary.phylog.cancor}}}

\examples{
data(SimulExample)
ex1.cancor <- cancor.phylog(SimulExample[,c(1,2,3,4,5)],SimulExample[,c(1,2,6,7,8)])
ex1.cancor
summary(ex1.cancor)
plot(ex1.cancor)
}
\keyword{multivariate}%-- one or more ...




