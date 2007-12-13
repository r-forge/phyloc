\name{check_phylo4d}
\alias{check_data}
\alias{attach_data}
\title{Validity checking for phylo4d objects}
\description{
  Basic checks on the validity of S4 phylogenetic plus data objects
}
\usage{
check_data(object,
		use.tip.names=TRUE,
		missing.tip.data=c("fail","OK","warn"),
		extra.tip.data=c("fail","OK","warn"),
		missing.tip.names=c("fail","OK","warn"),
		use.node.names=FALSE,
		missing.node.data=c("OK","warn","fail"),		
		extra.node.data=c("OK","warn","fail"),												
		missing.node.names=c("OK","warn","fail"))
}

attach_data(object,
    use.tip.names=TRUE,
		use.node.names=FALSE,
		use.edge.names=FALSE)
)

\arguments{
  \item{object}{A phylo4 object}
  \item{use.tip.names}{Use tip data names if present (default\=\code{TRUE})}
  \item{missing.tip.data}{Can tip data taxa be a subset of tree taxa\? (default\=\code{fail})}
  \item{extra.tip.data}{Can tip data taxa be a superset of tree taxa\? (default\=\code{fail})}
  \item{missing.tip.names}{Can tip data be supplied without taxa names\? (default\=\code{fail})}
  \item{use.node.names}{Use node data names if present (default\=\code{FALSE})}
  \item{missing.node.data}{Can node data taxa be a subset of tree taxa\? (default\=\code{fail})}
  \item{extra.node.data}{Can node data taxa be a superset of tree taxa\? (default\=\code{fail})}
  \item{missing.node.names}{Can node data be supplied without taxa names\? (default\=\code{fail})}
}
\value{
  \code{check_data} returns an error string (describing problems) or \code{TRUE} if everything is OK. If option \code{fail} is specified and a problem is encountered, returns an error string. If option \code{warn} is specified and a problem is encountered, returns \code{TRUE} plus a warning message.
  
  By default, tip data names are used to attach the tip data to the tree, must be present, and must match the set of tree taxa exactly. By default node data names are not used to attach the data to the tree, and thus node data names are not required.
  
  \code{attach_data} returns the phylo4 object with tipdata and nodedata (if present) attached to the phylo4 tree.
}
\details{
  \code{check_data} checks validity of a phylo4d object. For tipdata and nodedata, checks for presence of tip and node names, whether tip or node data taxa and tree taxa are an exact match, and whether tip or node data taxa are a superset or subset of the tree taxa.
  
  \code{attach_data} takes a phylo4d object that has been checked, and attaches the tipdata and nodedata to the object. This involves the data being subset to include only the taxa present in the tree, with missing data taxa assigned NA values. The data are sorted to match the order of the tip.labels of the tree to allow use with comparative analysis functions.
}
\author{Steven Kembel \email{skembel@berkeley.edu}}
\keyword{misc}
