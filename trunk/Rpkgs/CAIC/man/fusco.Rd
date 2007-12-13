\name{fusco.test}
\alias{fusco.test}
\alias{fusco.calc}
\alias{fusco.simulate}
\alias{fusco.null.tree}
\alias{summary.fusco}
\alias{print.fusco}
\alias{plot.fusco}


\title{Imbalance statistics using Fusco and Cronk's method.} 
\description{
INCOMPLETE - still in testing
}
\usage{
fusco.test(phy, rich = NULL, data = NULL, names = NULL, simulate=TRUE, sim.reps=1000, sim.CI = 0.95)
fusco.calc(phy, rich)
fusco.simulate <- function(obj, sim.reps=1000, sim.CI = 0.95)
fusco.null.tree <- function(nTips, nSpecies)
\method{print}{fusco}(obj, ...)
\method{summary}{fusco}(obj, ...)
\method{plot}{fusco}(x, smooth.correction=TRUE, breaks=10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{An object of class 'phylo'.}
  \item{data}{A data frame containing species richness values.}
  \item{names}{A variable in \code{data} identifying tip labels.}
  \item{rich}{A variable in \code{data} identifying species richness.}
  \item{simulate}{Use simulated Markov trees to generate confidence intervals on median I and quartile deviation.}
  \item{sim.reps}{Number of trees to simulate}
  \item{sim.CI}{Width of confidence intervals required.}
  \nTips{The number of tips in the tree to be simulated.}
  \nSpecies{The number of species distributed across those tips.}
  \item{obj}{An object of class 'fusco'.}
  \item{x}{A object of class 'fusco'.}
  \item{smooth.correction}{Apply the smoothing correction described in Appendix A of Fusco and Cronk (1995) to the histogram of nodal imbalance.}
  \item{breaks}{The number of bins to be used in the histogram of nodal imbalance.}

}
\details{
This is not complete and the results are not yet reconciled with the original application.

Informative nodes - bifurcating and with more than 3 species.
}
\value{
The function \code{fusco.test} produces an object of class 'fusco' containing:
	\item{observed}{A data frame of informative nodes showing nodal imbalance statistics.}
	\item{median}{The median value of I.}
	\item{qd}{The quartile deviation of I.}
	\item{tipsAsSpecies}{A logical indicating whether the tips of the trees were treated as species or higher taxa.}
	\item{nInformative}{The number of informative nodes.}
	\item{nSpecies}{The number of species distributed across the tips.}
	\item{nTips}{The number of tips.}
If \code{simulate} is TRUE, or the user calls \code{fusco.simulate} on a 'fusco' object, then the following are also present.
	\item{simulated}{A data frame of median and quartile deviation in I from each simulated tree.}
	\item{sim.median}{A vector of length 2 giving confidence intervals in median I from simulated trees.}
	\item{sim.qd}{A vector of length 2 giving confidence intervals in quartile deviation of I from simulated trees.}
	\item{sim.reps}{The number of replicates used.}
	\item{sim.CI}{The confidence levels used.}
}

\references{
Fusco, G. & Cronk, Q.C.B. (1995) A New Method for Evaluating the Shape of Large Phylogenies. J. theor. Biol. 175, 235-243
Purvis A.,  Katzourakis A. & Agapow, P-M (2002) Evaluating Phylogenetic Tree Shape: Two Modifications to Fusco & Cronk's Method. J. theor. Biol. 214, 93-103.
}
\author{David Orme, Andy Purvis}
\examples{
data(fusco)
fusco.test(fuscoBirdTree, fuscoBirdRich)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
