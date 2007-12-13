\name{brunch.contrast}
\alias{brunch.contrast}

\title{ Computes a BRUNCH contrast, if there is one, for a plum}
\description{
The taxa descended from a plum can be contrasted if they differ in the x-variable. This function tests whether 
they do and, if so, makes the contrast.}
\usage{
brunch.contrast(phylogeny, plum, arbitrary.bl = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ an object of class \code{"phylo"} }
  \item{plum}{ a plum; see \code{\link{find.plums}} }
  \item{arbitrary.bl}{ numeric; length of arbitrary branch used in treatment of polytomous plums. }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
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
\references{ Purvis, A. and Rambaut, A. 1995 \emph{Comput. Appls. Biosci.} \bold{11}, 247--251, for algorithms. }
\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{There is no error-checking}

\seealso{ \code{\link{find.plums}}, \code{\link{brunch}}, \code{\link{analyse.brunch.contrasts}}}

## The function is currently defined as
function(phylogeny, plum, arbitrary.bl=0){

#Computes a BRUNCH contrast if there is one.

branch.lengths<-rep(NA,length(plum$tips))
for (i in 1:length(plum$tips)) branch.lengths[i]<-get.branch.length(phylogeny, plum$tips[i])

if (var(plum[[4]])==0)
{
	#No variance in X within node, so compute average
	root.y<-sum((plum[[3]]/branch.lengths))/sum(1/branch.lengths)
	root.x<-sum((plum[[4]]/branch.lengths))/sum(1/branch.lengths)
	root.bl.adjustment<-1/sum(1/branch.lengths)
	y.contrast<-NA
	x.contrast<-NA
	variance<-NA
	to.drop<-plum$tips
  }

if (var(plum[[4]])!=0)
{
	#Variance in X within node, so compute contrast
	contrast.weights<-rep(NA,length(plum[[3]])) #construct vector to hold contrast weights
	mean.x<-mean(plum[[4]])
	contrast.weights<-ifelse((plum[[4]]>=mean(plum[[4]])),(1/sum(plum[[4]]>=mean.x)),(-1/sum(plum[[4]]<mean.x)))

	if (length(plum$tips)==2)
	{
	#Contrast at a bifurcation
	root.y<-NULL
	root.x<-NULL
	root.bl.adjustment<-NULL
	y.contrast<-sum(plum[[3]]*contrast.weights)
	x.contrast<-sum(plum[[4]]*contrast.weights)
	variance<-sum(branch.lengths)
	to.drop<-plum$tips
	}
	
	if (length(plum$tips)>2)
	{
	#Contrast at a polytomy
	root.y<-NULL
	root.x<-NULL
	root.bl.adjustment<-NULL
	va<-arbitrary.bl
	vb<-arbitrary.bl
	branch.lengths<-branch.lengths-arbitrary.bl
	a<-plum[[4]]>=mean(plum[[4]])
	ya<-sum((plum[[3]][a]/branch.lengths[a]))/sum(1/branch.lengths[a])
	xa<-sum((plum[[4]][a]/branch.lengths[a]))/sum(1/branch.lengths[a])
	yb<-sum((plum[[3]][!a]/branch.lengths[!a]))/sum(1/branch.lengths[!a])
	xb<-sum((plum[[4]][!a]/branch.lengths[!a]))/sum(1/branch.lengths[!a])
	va<-va+(1/sum(1/branch.lengths[a]))
	vb<-vb+(1/sum(1/branch.lengths[!a]))
	
	y.contrast<-ya-yb
	x.contrast<-xa-xb
	variance<-va+vb
	to.drop<-plum$tips
	}
  }

contrast<-list(nodal.y=root.y, nodal.x=root.x, bl.adjustment=root.bl.adjustment, y.contrast=y.contrast, x.contrast=x.contrast,
	variance=variance, to.drop=to.drop)
	
contrast
  }
}
\keyword{ htest}
\keyword{ regression}