\name{as.phylo4.phylo}
\alias{as.phylo4.phylo}
\alias{as.phylo.phylo4}
\title{Translate between phylo and phylo4 classes}
\description{
  Translates between phylo (S3, ape) and phylo4 (S4) classes of
  phylogenetic trees
}
\usage{
as.phylo4.phylo(x)
as.phylo.phylo4(x)
}
\arguments{
  \item{x}{A phylogenetic tree of the appropriate class}
}
\value{
  A phylogenetic tree object of the appropriate class
}
\author{Ben Bolker}
\examples{
    library(ape)
  example(read.tree)
  p1 = tree.owls
  P1 = as.phylo4.phylo(tree.owls)
  P1
  identical(p1,as.phylo.phylo4(P1))
  p1$edge ## works!

}
\keyword{misc}

