\name{phylo.overlap}
\alias{phylo.overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the phylogenetic diversity overlap between two assemblages.}
\description{
Calculates phylogenetic diversity (PD) as Total Branch Length for either the species shared by two assemblages or the species unique to the first assemblage.
}
\usage{
phylo.overlap(assemblage.one, assemblage.two, cm, method = c("shared", "unique"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assemblage.one}{A vector identifying the tips present in the first assemblage.}
  \item{assemblage.two}{A vector identifying the tips present in the first assemblage.  }
  \item{cm}{ An object of class 'clade.matrix'}
  \item{method}{One of either 'shared' or 'unique', giving the PD that is either shared between the two assemblages or is unique to the first.}
}
\details{
Note that phylo.overlap(1, 2, cm, "shared") is equivalent to phylo.overlap(2, 1, cm, "shared") but that phylo.overlap(1, 2, cm, "unique") is not the same as phylo.overlap(2, 1, cm, "unique").
}
\value{
A numeric vector of length one giving the PD.
}
\author{David Orme }

\seealso{ \code{\link{pd.calc}}}
\examples{
}
\keyword{ utilities }% at least one, from doc/KEYWORDS
