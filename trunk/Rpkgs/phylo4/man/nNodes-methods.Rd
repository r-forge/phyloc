\name{phylo4-accessors}
\docType{methods}
\alias{nNodes}
\alias{nNodes-methods}
\alias{nNodes,phylo4-method}
\alias{nTips}
\alias{nTips-methods}
\alias{nTips,phylo4-method}
\alias{nTips,ANY-method}
\alias{labels,phylo4-method}
\alias{edges}
\alias{edges-methods}
\alias{edges,phylo4-method}
\alias{nEdges}
\alias{nEdges-methods}
\alias{nEdges,phylo4-method}
\alias{hasEdgeLength}
\alias{hasEdgeLength-methods}
\alias{hasEdgeLength,phylo4-method}
\alias{EdgeLength}
\alias{EdgeLength-methods}
\alias{EdgeLength,phylo4-method}
\alias{isRooted}
\alias{isRooted-methods}
\alias{isRooted,phylo4-method}
\alias{RootEdge}
\alias{RootEdge-methods}
\alias{RootEdge,phylo4-method}
\alias{hasNodeLabels}
\alias{hasNodeLabels-methods}
\alias{hasNodeLabels,phylo4-method}
\alias{NodeLabels}
\alias{NodeLabels-methods}
\alias{NodeLabels,phylo4-method}
\alias{tdata}
\alias{tdata-methods}
\alias{tdata,phylo4d-method}
\title{Methods for S4 phylog classes}
\description{
   Generic methods for phylogenetic trees represented as S4 classes
}
\section{Methods}{
\describe{
  \item{nTips}{\code{signature(object="phylo4")}: number of tips}
  \item{labels}{\code{signature(object = "phylo4")}: tip labels}
  \item{nNodes}{\code{signature(object="phylo4")}: number of internal nodes}
  \item{hasNodeLabels}{\code{signature(object = "phylo4")}: whether
    tree has (internal) node labels}
  \item{NodeLabels}{\code{signature(object = "phylo4")}: internal
    node labels}
  \item{nEdges}{\code{signature(object = "phylo4")}: number of edges}
  \item{edges}{\code{signature(object = "phylo4")}: edge matrix}
  \item{hasEdgeLength}{\code{signature(object = "phylo4")}: whether
    tree has edge (branch) lengths}
  \item{EdgeLength}{\code{signature(object = "phylo4")}: edge (branch)
    lengths (or NULL if none)}
  \item{isRooted}{\code{signature(object = "phylo4")}: whether tree is
    rooted (i.e. has explicit root edge defined \emph{or} root node has
    <= 2 descendants)}
  \item{RootEdge}{\code{signature(object = "phylo4")}: root edge}
  \item{tdata}{\code{signature(object = "phylo4d")}: trait data}
}
}
\section{usage}{
  \S4method{tdata}{phylo4d}(x,which="tip",\dots)
}
\arguments{
  \item{which}{Which data to extract: "tip" (tips only), "node"
    (internal nodes only), "allnode" (tips and internal nodes), "edge"
    (edges) }
}
\keyword{methods}
