\name{root.to.tip}
\alias{root.to.tip}
\title{Compute distance between the root of a phylogeny and a specified tip}
\description{
  Sums the lengths of branches between the root of a phylogeny and a specified tip.
}
\usage{
root.to.tip(phylogeny, tip.name)
}
\arguments{
  \item{phylogeny}{ An object of class phylo }
  \item{tip.name}{ An element of the phylogeny's tip.label }
}
\details{
  Called by \code{get.root.to.tip}.}
\value{
  numeric
}
\section{Warning}{

  \code{root.to.tip} does not yet check that it has been passed valid arguments.

}

\author{ Andy Purvis }

\seealso{ \code{\link{get.root.to.tip}}}

## The function is currently defined as
function(phylogeny, tip.name){
#Compute a single root-to-tip value for taxon tip.name
n.label<-which(phylogeny$tip.label==tip.name)
n.edge<-which(phylogeny$edge[,2]==as.character(n.label))
sum.so.far<-phylogeny$changes[n.edge]
while(phylogeny$edge[n.edge,1]!="-1")
{
	n.edge<-which(phylogeny$edge[,2]==phylogeny$edge[n.edge,1])
	sum.so.far<-sum.so.far+phylogeny$changes[n.edge]
  }
return(sum.so.far)
  }
}
\keyword{ manip }
