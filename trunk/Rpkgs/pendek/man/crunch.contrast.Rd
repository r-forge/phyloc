\name{crunch.contrast}
\alias{crunch.contrast}
\title{ Computes a CRUNCH contrast; NOT YET COMPLETED }
\description{
This function will compute a CRUNCH-style contrast at each node in a phylogeny, without having 
recourse to \code{pic}. It is not yet complete.}
\usage{
crunch.contrast(phylogeny, plum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{"phylo"} }
  \item{plum}{ A plum, from \code{\link{find.plums}} }
}

}
\value{
A list with the following elements:
  \item{nodal.y }{weighted average of y values if no contrast computed, otherwise NA}
  \item{nodal.x }{weighted average of x values if no contrast computed, otherwise NA}
  \item{bl.adjustment }{adjustment required to parental branch length if no contrast computed, otherwise NA}
  \item{y.contrast }{contrast in y-variable if one is computed, otherwise NA}
  \item{x.contrast }{contrast in x-variable if one is computed, otherwise NA}
  \item{variance }{expected variance of contrast if one is computed, otherwise NA}
  \item{to.drop }{taxa used in the comparison}
}
\references{ Pagel, M. D. 1992 \emph{J. theor. Biol.} \bold{156}, 431--442 for algorithms. }
\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{There is no error-checking}

\seealso{ \code{\link{find.plums}}, \code{\link{brunch.contrast}}}


## The function is currently defined as
function(phylogeny, plum){

#Computes a CRUNCH contrast.

branch.lengths<-rep(NA,length(plum$tips))
for (i in 1:length(plum$tips)) branch.lengths[i]<-get.branch.length(phylogeny, plum$tips[i])


contrast.weights<-rep(NA,length(plum[[3]])) #construct vector to hold contrast weights
mean.x<-mean(plum[[4]])
contrast.weights<-ifelse((plum[[4]]>=mean(plum[[4]])),(1/sum(plum[[4]]>=mean.x)),(-1/sum(plum[[4]]<mean.x)))

if (length(plum$tips)==2)
{
	#Contrast at a bifurcation
	root.y<-sum((plum[[3]]/branch.lengths))/sum(1/branch.lengths)
	root.x<-sum((plum[[4]]/branch.lengths))/sum(1/branch.lengths)
	root.bl.adjustment<-1/sum(1/branch.lengths)
	y.contrast<-sum(plum[[3]]*contrast.weights)
	x.contrast<-sum(plum[[4]]*contrast.weights)
	variance<-sum(branch.lengths)
	to.drop<-plum$tips
  }
	
if (length(plum$tips)>2)
	{
	#Contrast at a polytomy
	#At present, does nothing different from bifurcations - this needs more work!!!
	root.y<-sum((plum[[3]]/branch.lengths))/sum(1/branch.lengths)
	root.x<-sum((plum[[4]]/branch.lengths))/sum(1/branch.lengths)
	root.bl.adjustment<-1/sum(1/branch.lengths)
	y.contrast<-sum(plum[[3]]*contrast.weights)
	x.contrast<-sum(plum[[4]]*contrast.weights)
	variance<-sum(branch.lengths)
	to.drop<-plum$tips
	}


contrast<-list(nodal.y=root.y, nodal.x=root.x, bl.adjustment=root.bl.adjustment, y.contrast=y.contrast, x.contrast=x.contrast,
	variance=variance, to.drop=to.drop)
	
contrast
  }
}
\keyword{regression}