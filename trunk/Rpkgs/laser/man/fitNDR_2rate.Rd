\name{fitNDR_2rate.Rd}
\alias{fitNDR_2rate}

\title{ Fit 2-rate diversification model to combined phylogenetic/taxonomic data  }
\description{
  Finds likelihoods and parameter estimates for combined phylogenetic/taxonomic data under a two-rate diversification model. 
}    
\usage{
fitNDR_2rate(phy, eps=0, combined=TRUE, rate.decrease=FALSE)

}


\arguments{
  \item{phy}{ a class 'phylo' phylogenetic tree with an additional component 'phenotype'; \code{phy$phenotype} specifies the number of tips per terminal taxon.  The phenotype data is added with \code{getTipData}}
  \item{eps}{ the relative extinction rate, or \code{mu / lambda} }
  \item{ combined }{ Should likelihoods be calculated using both combined phylogenetic and taxonomic data, or taxonomic data only? }
  \item{rate.decrease}{Constrain the analysis such that the highest diversification rate occurs in the tree bipartition containing the root node?}
}
\details{
  \code{fitNDR_2rate} fits a 2-speciation rate model to combined phylogenetic/taxonomic data. Suppose you have a higher-level phylogeny for some group of organisms (e.g., beetle families), where you also know the approximate species diversities for each terminal. \code{fitNDR_2rate} assumes that at some point in the tree, an ancestral speciation rate \code{lambda1} shifts to a new speciation rate \code{lambda2}. 
  
  The model is fitted by iteratively splitting the tree at each node and fitting a birth-death model to each of the resulting bipartitions. Thus, for each node, you obtain (i) the likelihood of a rate shift at that position in the tree, and (ii) the estimated speciation rates for each bipartition. The function returns a dataframe giving the likelihoods of rate shift at each node as well as the parameter estimates. 
  
    \code{eps} allows you to estimate the speciation rate under any assumed (constant) relative extinction rate, where the relative extinction rate is \code{mu / lambda }. 
  
  Specifying \code{combined = TRUE} estimates \code{lambda} from combined taxonomic and phylogenetic data. Thus, the likelihood is a function of the internal branching structure of the tree and the species richness/taxonomic data. If \code{combined=FALSE}, likelihoods are estimated from the species richness/taxonomic data only, and the internal branching structure of the tree does not contribute.
  
  Option \code{rate.decrease} fits a model where the highest speciation rate must occur in the tree bipartition containing the root node.

}
\value{
  a dataframe with the following components:
  \item{node}{ The node defining the tree bipartition }
  \item{LH }{ The log-likelihood at the maximum}
  \item{aic }{the Akaike Information Criterion}
  \item{r.1}{ the net diversification rate,\code{lambda - mu}, for the partition containing the root node}
  \item{lambda.1}{the speciation rate for the for the partition containing the root node}
  \item{LH.1}{the log-likelihood for partition 1 at the maximum}
   \item{r.2}{ the net diversification rate,\code{lambda - mu}, for the partition NOT containing the root node, e.g., partition 2}
  \item{lambda.2}{the speciation rate for the for the partition 2} 
  \item{LH.2}{the log-likelihood for partition 2 at the maximum}
  \item{eps }{ the value of \code{eps} used}
  
}
\references{ 

Rabosky, D. L., S. C. Donnellan, A. L. Talaba, and I. J. Lovette. 2007. Exceptional among-lineage variation in diversification rates during the radiation of Australia’s largest vertebrate clade. Proc. Roy. Soc. Lond. Ser. B 274:2915-2923.


}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{fitNDR_1rate}}, \code{\link{getTipdata}}, \code{\link{lambda.stem.ml}}  }
\examples{

}
\keyword{ models }
