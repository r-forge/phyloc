\name{get.branch.length}
\alias{get.branch.length}

\title{ Find branch length leading to a specified phylogenetic tip }
\description{
Returns the branch length leading only to the specified phylogenetic tip}
\usage{
get.branch.length(phylogeny, tip.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ an object of class \code{"phylo"} }
  \item{tip.name}{ character string; an element of \code{phylogeny$tip.label}}
}

}
\value{
numeric}
\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking}


## The function is currently defined as
function(phylogeny, tip.name){

#Returns the edge length leading to a particular tip in the phylogeny

edge.label<-which(phylogeny$tip.label==tip.name)
edge.index<-which(phylogeny$edge[,2]==as.character(edge.label))

return(phylogeny$edge.length[edge.index])
  }
}
\keyword{ manip}