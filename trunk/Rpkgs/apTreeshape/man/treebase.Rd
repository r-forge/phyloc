\name{treebase}
\alias{primates}
\alias{carnivora.treeshape}

\title{Phylogenetic trees from the TreeBASE database}

\description{
This set contains two trees extracted from the TreeBASE Database.\cr
\code{carnivora.treeshape}: \cr
This data set describes the phylogenetic relationships of carnivora as reported by Bininda-Emonds OR, Gittleman JL, Purvis A.. They inferred this phylogeny combining phylogenetic information in a supertree construction.\cr  
\code{primates}:\cr
This data set describes the phylogenetic relationships of the primates as reported by Hayasaka K., Gojobori T. and Horai S.. They inferred this phylogeny from the evolution of mitochondrial DNA.
  
}

\usage{
data(carnivora.treeshape)
data(primates)
}

\format{
These data are stored as an object of class \code{"treeshape"}.
}

\source{
Bininda-Emonds, O. R., Gittleman, J. L. and Purvis, A. (1999) Building large trees by combining phylogenetic information: a complete phylogeny of the extant Carnivora (Mammalia). \emph{Biological Reviews}, \bold{74}, 143 -- 175. for more details about \code{"carnivora.treeshape"}.
  
Hayasaka, K., Gojobori, T. and Horai, S. (1988) Molecular phylogeny and evolution of primate mitochondrial DNA. \emph{Molecular Biology and Evolution}, \bold{5}, 626 -- 644, for more details about \code{"primates"}. 
}

\examples{
  
# Example tree in Treeshape format
data(carnivora.treeshape) 
plot(carnivora.treeshape)

data(primates)
plot(primates)
}

\keyword{datasets}
