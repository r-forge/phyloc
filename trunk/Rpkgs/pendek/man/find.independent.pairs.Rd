\name{find.independent.pairs}
\alias{find.independent.pairs}

\title{Decompose a fully binary phylogeny into independent pairs of species }
\description{
This function decomposes a fully binary phylogeny into independent pairs of species, using a 
version of the BRUNCH algorithm simplified to assume that the phylogeny is fully binary.}
\usage{
find.independent.pairs(phylogeny)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ An object of class \code{"phylo"}, assumed to be fully binary. }
}
\details{
Proceeds from root to tips, finding cherries and pruning them from the tree until none remain.}
\value{
An n x 2 matrix of terminal names taken from \code{phylogeny$tip.label}.}


\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{There is no error-checking}

\seealso{ \code{\link{find.cherries}}}


## The function is currently defined as
function(phylogeny){

#Finds all n phylogenetically independent pairs of taxa in a fully binary phylogeny
#Returns an n x 2 matrix of the tip names

n<-floor(length(phylogeny$tip.label)/2) #Number of comparisons
pairs<-matrix(NA, nrow=n ,ncol=2)
count<-0
while (count<n)
{
	some.pairs<-find.cherries(phylogeny)
	how.many<-length(some.pairs[,1])
	pairs[(count+1):(count+how.many),]<-some.pairs
	count<-count+how.many
	to.drop<-some.pairs
	if (length(phylogeny$tip.label)-1!=length(to.drop)) phylogeny<-drop.tip(phylogeny, to.drop)
  }

pairs
  }
}
\keyword{ manip}
\keyword{ htest}