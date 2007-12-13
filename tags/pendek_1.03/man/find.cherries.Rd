\name{find.cherries}
\alias{find.cherries}
\title{Finds all cherries (pairs of sister terminals) in a phylogeny}
\description{
Cherries are pairs of terminal taxa in a phylogeny that are each other's closest relatives, 
at least among the terminal taxa in the phylogeny.  This function finds all such pairs.  It is
used by various other functions in the package.}
\usage{
find.cherries(phylogeny)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{"phylo"}}
}
}
\value{
An n x 2 matrix, listing the \code{tip.label}s of terminals within each cherry.}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking on the phylogeny passed to the function}

\seealso{\code{\link{find.plums}}, \code{\link{brunch.contrast}} }

## The function is currently defined as
function(phylogeny){

#Finds the n cherries in a fully binary phylogeny
#Returns an n x 2 matrix with the tip names

ancestral.nodes<-which(table(phylogeny$edge[, 1][phylogeny$edge[, 2] > 0])==2)
ancestors<-names(ancestral.nodes)
pairs<-matrix(NA,nrow=length(ancestral.nodes),ncol=2)

for (i in 1:length(ancestors))
{
	sisters<-phylogeny$edge[,2][which(phylogeny$edge[,1]==ancestors[i])]
	pairs[i,1]<-phylogeny$tip.label[as.numeric(sisters[1])]
	pairs[i,2]<-phylogeny$tip.label[as.numeric(sisters[2])]
  }
pairs
  }
}
\keyword{manip}