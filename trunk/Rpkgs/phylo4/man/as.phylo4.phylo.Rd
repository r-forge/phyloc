\name{as.phylo4.phylo}
\alias{as.phylo4.phylo}
\alias{as.phylo.phylo4}
\alias{as.multi.tree.phylo4}
\alias{as.phylo4.multi.tree}
\alias{as.phylo4}
\title{Translate between phylo and phylo4 classes}
\description{
  Translates between phylo (S3, ape) and phylo4 (S4) classes of
  phylogenetic trees
}
\usage{
as.phylo4(x,\dots)
\method{as.phylo4}{phylo}(x,\dots)
\method{as.phylo4}{multi.tree}(x,\dots)
as.phylo.phylo4(x)
as.multi.tree.phylo4(x)
}
\arguments{
  \item{x}{A phylogenetic tree of the appropriate class}
  \item{\dots}{additional arguments, for future generic-compatibility}
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
}
\keyword{misc}

