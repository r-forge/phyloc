\name{pandit}
\alias{cytochromc}
\alias{rhodopsin}

\title{Phylogenetic trees from the Pandit database}

\description{
This set contains two trees extracted from the Pandit Database,\code{cytochromc} and \code{rhodopsin}.
}

\usage{
data(cytochromc)
data(rhodopsin)
}

\details{ 
Cytochromes c (cytC) can be defined as electron-transfer proteins.This data set describes the phylogenetic relationships of cytC of mitochondria and bacteria.\cr 
G-protein-coupled receptors, GPCRs, constitute a vast protein family that encompasses a wide range of functions (including various autocrine, paracrine and endocrine processes). The rhodopsin-like GPCRs themselves represent a widespread protein family that includes hormone, neurotransmitter and light receptors.
}
 
\source{
Ambler, R. P. (May, 1991) Sequence variability in bacterial cytochromes c. \emph{Biochimica et Biophysica Acta}, \bold{1058}, 42-7, for more details about bacterial cytochromes c. 
  
Attwood T. K. and Findlay J. B. (Feb., 1993) Design of a discriminating fingerprint for G-protein-coupled receptors. \emph{Protein Engineering}, \bold{6}, 176-67, for more details about G-protein-coupled receptors.
}

\format{
These data are stored as objects of class \code{"treeshape"}.
}

\examples{

## Example tree in Treeshape format
data(cytochromc) 
plot(cytochromc)
summary(cytochromc) 

data(rhodopsin) 
plot(rhodopsin)
}

\keyword{datasets}
