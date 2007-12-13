\name{phylo4d}
\alias{phylo4d}
\title{Combine tree with data}
\description{
  Merges a phylogenetic tree (class \code{phylo4}) with
  a data frame to create a combined object of class \code{phylo4d}
}
\usage{
phylo4d(tree, data, which = "tip")
}
\arguments{
  \item{tree}{an object of class \code{phylo4}}
  \item{data}{a data frame for node data (either tips, internal nodes,
    or all)}
%  \item{edgedata}{a data frame for edge traits}
  \item{which}{character: "tip", "node", or "all" depending on
    which data are specified}
}
\value{
  An object of class \code{phylo4d}
}
\author{Ben Bolker}
\note{
  Checking on matches will be done by the validity checker
  (label matches between data and tree tips, number of rows
  of data vs. number of nodes/tips/etc.)
}
\examples{
library(ape)
example(read.tree)
phylo4d(as(tree.owls.bis,"phylo4"),data=data.frame(wing=1:3))
}
\keyword{misc}
