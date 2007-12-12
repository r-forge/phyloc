\name{get.root.to.tip}
\alias{get.root.to.tip}
\title{Computes the distance between a phylogeny's root and each of its tips.}
\description{
  For each tip in the phylogeny, computes the sum of the branch lengths between the tip
  and the root. The root-to-tip distances are returned as a vector in the same order as
  as the tips.
}
\usage{
get.root.to.tip(phylogeny)
}
\arguments{
  \item{phylogeny}{ An object of class phylo }
}

\details{
  Calls \code{root.to.tip} for each of the tips in turn.
}
\value{
  A numeric vector, in the same order as the tip.labels of the phylogeny.
}
\author{ Andy Purvis }

\section{Warning}{

  \code{get.root.to.tip} currently does not check whether it has been passed a valid phylo object.


}

\seealso{\code{\link{root.to.tip}}}

## The function is currently defined as
function(phylogeny){
#Compute root to tip sum of character changes for each tip
r2t.changes<-rep(NA,length(phylogeny$tip.label))
for(i in 1:length(phylogeny$tip.label))
{	
	r2t.changes[i]<-root.to.tip(phylogeny,phylogeny$tip.label[i])
  }
return(r2t.changes)
  }
}
\keyword{ manip }