\name{mlphylo}
\alias{mlphylo}
\alias{logLik.phylo}
\alias{deviance.phylo}
\alias{AIC.phylo}
\title{Estimating Phylogenies by Maximum Likelihood}
\usage{
mlphylo(model = DNAmodel(), x, phy, search.tree = FALSE, quiet = FALSE)
\method{logLik}{phylo}(object, ...)
\method{deviance}{phylo}(object, ...)
\method{AIC}{phylo}(object, ..., k = 2)
}
\arguments{
  \item{model}{an object of class \code{"DNAmodel"} giving the model to
    be fitted.}
  \item{x}{a matrix, a data frame, or a list giving the (aligned) DNA
    sequence data.}
  \item{phy}{an object of class \code{"phylo"} giving the (initial)
    tree.}
  \item{search.tree}{a logical specifying whether to search for the best
    tree (defaults to FALSE) (not functional for the moment).}
  \item{quiet}{a logical specifying whether to display the progress of
    the analysis.}
  \item{object}{an object of class \code{"phylo"}.}
  \item{k}{a numeric value giving the penalty per estimated parameter;
    the default is \code{k = 2} which is the classical Akaike
    information criterion.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  \code{mlphylo} estimates a phylogenetic tree by maximum likelihood
  given a set of DNA sequences. The model of evolution is specified with
  the function \code{\link{DNAmodel}}.

  \code{logLik}, \code{deviance}, and \code{AIC} are generic functions
  used to extract the log-likelihood, the deviance (-2*logLik), or the
  Akaike information criterion of a tree. If no such values are
  available, \code{NULL} is returned.
}
\details{
  The present version is a pre-alpha release. All comments, suggestions,
  bug reports, are warmly welcome.

  The model specified by \code{\link{DNAmodel}} is fitted using the
  standard ``pruning'' algorithm of Felsenstein (1981). An algorithm for
  the estimation of tree topology is under development, and will be
  released when ready.

  The implementation of the inter-sites variation in substitution rates
  follows the methodology developed by Yang (1994).

  The difference among partitions is parametrized with a contrast
  parameter (denoted \eqn{\xi}{xi}) that specifies the contrast in mean
  susbtitution rate among the partitions. This methodology is inspired
  from one introduced by Yang (1996).

  The substitution rates are indexed column-wise in the rate matrix: the
  first rate is set to one.
}
\note{
  For the moment, it is not possible to estimate neither branch lengths,
  nor the topology with \code{mlphylo}: this is still in development and
  will be released when ready. The function may estimate all other
  parameters: substitution rates, shape (\eqn{\alpha}{alpha}) of the
  inter-sites variation in substitution rates, the proportion of
  invariants, and the ``contrast'' parameter (\eqn{\xi}{xi}) among
  partitions.

  Alternative topologies can also be compared using likelihood-ratio
  tests (LRTs) or AICs.
}
\value{
  an object of class \code{"phylo"} with branch lengths as estimated by
  the function. There are two additional attributes:

  \item{loglik}{the maximum log-likelihood.}
  \item{para}{the estimated parameters for each partition.}
}
\references{
  Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a
  maximum likelihood approach. \emph{Journal of Molecular Evolution},
  \bold{17}, 368--376.

  Yang, Z. (1994) Maximum likelihood phylogenetic estimation from DNA
  sequences with variable rates over sites: approximate methods.
  \emph{Journal of Molecular Evolution}, \bold{39}, 306--314.

  Yang, Z. (1996) Maximum-likelihood models for combined analyses of
  multiple sequence data. \emph{Journal of Molecular Evolution},
  \bold{42}, 587--596.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{DNAmodel}}, \code{\link{nj}}, \code{\link{read.dna}},
  \code{\link{summary.phylo}}
}
\keyword{models}
