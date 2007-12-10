\name{prcomp.phylog}
\alias{prcomp.phylog}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Principal Components Analysis from Simulated Data Sets}
\description{
 Performs a principal component analyses on a set of simulated data and
 return the eigenvalues.
}
\usage{
\method{prcomp}{phylog}(x, max.num=0, exclude.tips=NULL, lapply.size=100,
              center=TRUE, scale=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ the columns from a data set returned from read.sim.data
   that you want to use in the PCA. The first column MUST be sim.counter
 and the second Tips.}
 \item{max.num}{ if different from 0, the maximum number of simulations
   to analyze. }
 \item{exclude.tips}{ an optional vector giving the names of tips to
   exclude from the analyses.}
 \item{lapply.size}{ a tuning parameter that can affect the speed of
   calculations; see Details in phylog.lm.}
 \item{center}{ should the data be centered before analyses? defaults to
 yes; see prcomp.}
 \item{scale}{ should the data be scaled before the analyses? defaults
   to yes; see prcomp. }
 \item{...}{Not used.}
}
\value{ A list (of class phylog.prcomp) with components
  \item{call}{the call to the function}
  \item{Eigenvalues}{all the eigenvalues from the PCA.  The one with
    sim.counter=0 corresponds to the original (''real'') data.}
}

\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
  
  Krzanowski, W. J. (1990) \emph{Principles of multivariate analysis}
  Oxford University Press.
  
  Morrison, D. F. (1990) \emph{Multivariate statistcal methods, 3rd ed.} McGraw-Hill.
}
\author{Ramon Diaz-Uriarte and Theodore Garland, Jr.}



\seealso{ \code{\link{read.sim.data}}, \code{\link{summary.phylog.prcomp}}}

\examples{
data(SimulExample)
ex1.prcomp <- prcomp.phylog(SimulExample[,-11]) # 11th col. is a factor
options(digits=5)
ex1.prcomp
summary(ex1.prcomp)
plot(ex1.prcomp)
}
\keyword{multivariate}%-- one or more ...
