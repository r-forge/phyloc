\name{fitNDR_1rate.Rd}
\alias{fitNDR_1rate}

\title{ Estimate diversification rate from combined taxonomic/phylogenetic data  }
\description{
  Find maximum likelihood estimate of net diversification rate from combined taxonomic/phylogenetic data.
}    
\usage{
fitNDR_1rate(phy, eps=0, combined=TRUE)

}


\arguments{
  \item{phy}{ a class 'phylo' phylogenetic tree with an additional component 'phenotype'; \code{phy$phenotype} specifies the number of tips per terminal taxon.  The phenotype data is added with \code{getTipdata}}
  \item{eps}{ the relative extinction rate, or \code{mu / lambda} }
  \item{ combined }{ Should likelihoods be calculated using both combined phylogenetic and taxonomic data, or taxonomic data only? }
}
\details{
  \code{fitNDR_1rate} finds the maximum likelihood estimate of the net diversification rate using phylogenetic and taxonomic data. This method is best applied to phylogenetic trees for higher taxonomic levels where you have incomplete sampling but know the approximate species diversities of each terminal taxon in the tree (e.g., a phylogenetic tree of arthropod families).
  
  \code{eps} allows you to estimate the speciation rate under any assumed (constant) relative extinction rate, where the relative extinction rate is \code{mu / lambda }. 
  
  The \code{combined} argument asks whether you wish to use an estimator based on (i) combined taxonomic and phylogenetic data, or (ii) taxonomic data only. There are good reasons for trying both approaches: when \code{combined = TRUE}, both the taxonomic/species richness and phylogenetic backbone of a tree contribute heavily to the overall likelihood estimate.  This occurs in spite of the fact that most of the species diversity is actually in the tips. When \code{combined = TRUE}, the internal phylogenetic structure of the tree does not contribute to the likelihood; you are finding the maximum likelihood estimate of the speciation rate \code{lambda} from the species richness data plus stem clade ages alone. I recommend checking estimates under both \code{combined =TRUE} and \code{combined=FALSE}

}
\value{
  a dataframe with the following components:
  \item{LH }{ The log-likelihood at the maximum}
  \item{aic }{the Akaike Information Criterion}
  \item{r }{ the net diversification rate,\code{lambda - mu}, at the maximum}
  \item{lambda }{ the ML estimate of the speciation rate}
  \item{eps }{ the value of \code{eps} used}
  
}
\references{ 

Rabosky, D. L., S. C. Donnellan, A. L. Talaba, and I. J. Lovette. 2007. Exceptional among-lineage variation in diversification rates during the radiation of Australia’s largest vertebrate clade. Proc. Roy. Soc. Lond. Ser. B 274:2915-2923.


}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{fitNDR_2rate}}, \code{\link{getTipdata}}, \code{\link{lambda.stem.ml}} }
\examples{

}
\keyword{ models }
