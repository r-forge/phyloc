\name{phylo4-class}
\docType{class}
\alias{phylo4-class}
\alias{phylo4d-class}
\alias{$,phylo4-method}
\alias{print,phylo4-method}
\alias{show,phylo4-method}
\alias{summary,phylo4-method}
\title{phylo4 and extended classes}
\description{Classes for phylogenetic trees}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("phylo4", ...)
  Translation from other phylogenetic packages ...}.
}
\section{Slots}{
  \describe{
    \item{\code{edge}:}{Matrix of edges}
    \item{\code{edge.length}:}{Numeric vector of edge lengths}
    \item{\code{Nnode}:}{Integer: number of nodes}
    \item{\code{tip.label}:}{Character vector of tip labels}
    \item{\code{root.edge}:}{Integer: root edge (NA if none)}
  }
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "phylo4")}: extract a slot}
    \item{print}{\code{signature(x = "phylo4")}: print method}
    \item{show}{\code{signature(object = "phylo4")}: show method }
    \item{summary}{\code{signature(object = "phylo4")}: summary method}
    }
  }
  \section{usage}{
    \S4method{print}{phylo4}(x,printlen,...)
  }
  \arguments{
    \item{printlen}{Maximum number of tip labels to print}
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
\keyword{classes}
