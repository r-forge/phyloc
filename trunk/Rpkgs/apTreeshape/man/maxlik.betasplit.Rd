\name{maxlik.betasplit}
\alias{maxlik.betasplit}


\title{ Maximum likelihood of beta in the Beta-splitting model }

\description{ The function finds the beta value that maximizes the likelihood in the Beta-splitting model. Beta=0 corresponds
to Yule trees, Beta<0 corresponds to trees more unbalanced than Yule trees and Beta>0 corresponds to trees more balanced than Yule trees. Confidence intervals may also be provided.
}

\usage{
maxlik.betasplit(phylo, up = 10, remove.outgroup = FALSE, confidence.interval = "none", prob.conf.inter = 0.95, size.bootstrap = 100)
}


\arguments{
  \item{phylo}{An object of class '"treeshape"' or '"phylo"' on which the likelihood is computed.}
  \item{up}{ numeric. The Beta value that maximizes the likelihood is searched between -2 and \code{up} }
  \item{remove.outgroup}{ logical. Should one or two outgroups be removed before computing the likelihood }
  \item{confidence.interval}{ The method to be used for computing confidence intervals. See *Details*. }
  \item{prob.conf.inter}{ numeric. Probability associated with the confidence interval }
  \item{size.bootstrap}{ number. Number of bootsrap replicates }
}
\details{
The beta-splitting model has been introduced by Aldous to simulate trees with different tree balance.

Beta=0 corresponds to Yule trees.

Beta<0 corresponds to trees more unbalanced than Yule trees, Beta=-3/2 corresponds to the PDA model.

Beta>0 corresponds to trees more balanced than Yule trees.

By default, \code{confidence.interval}='"none"' and no confidence interval is computed.

When \code{confidence.interval}=\code{bootstrap_param}, parametric bootstraping is performed to compute
a confidence interval. Shall be used when the number of tips is small (<50).

When \code{confidence.interval}=\code{chi_square}, a chi-square approximatiom is used for
computing a confidence interval. Shall be used when the number of tips is large (>50).

}

\value{
  \item{max_lik}{The Beta value that maximizes the likelihood}
  \item{conf_interval}{A confidence interval for \code{max_lik}}
}

\references{
Aldous, D. J. (1996) Probability Distributions on Cladograms pp.1-18 of Random Discrete Structures eds D. Aldous and R. Pemantle, IMA Volumes Math. Appl. 76. 

Aldous, D. J. (2001) Stochastic Models and Descriptive Statistics for Phylogenetic Trees, from Yule to Today. Statistical Science, *16*, 23 - 24.

Blum, M.G.B. and Francois, O. Which random processes describe the Tree of Life? A large-scale study of phylogenetic tree imbalance. Systematic Biology *55*, 685-691, 2006.
}

\author{ michael.blum@imag.fr }

\seealso{ 
\code{\link{sackin}},\code{\link{sackin.test}},\code{\link{colless}},\code{link{colless.test}} 
}
\examples{
library(ape)
library(Bhat)
tree.pda<-rpda(50)
maxlik.betasplit(tree.pda,confidence.interval="none")
maxlik.betasplit(tree.pda,confidence.interval="chi_square")
maxlik.betasplit(tree.pda,confidence.interval="bootstrap_param")
}

\keyword{htest}