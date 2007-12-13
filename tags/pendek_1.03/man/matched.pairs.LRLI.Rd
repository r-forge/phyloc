\name{matched.pairs.LRLI}
\alias{matched.pairs.LRLI}

\title{ LRLI differences on phylogenetically independent matched pairs }
\description{
Computes and returns log(Rate) vs log(Interval) (LRLI) for phylogenetically independent 
matched pairs}
\usage{
matched.pairs.LRLI(phylogeny, phylogeny.matrix, pairs, data.matrix, plot = "n")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{phylo}}
  \item{phylogeny.matrix}{ A matrix of phylogenetic distances between pairs of terminals, as 
  constructed by \code{dist.phylo}. }
  \item{pairs}{ n x 2 matrix of tip labels, as produced by \code{find.independent.pairs}}
  \item{data.matrix}{ matrix of data differences, as constructed by \code{as.matrix(dist())}. }
  \item{plot}{ is a plot required? (Default = n) In fact, typing a "y" won't help anyway }
}
\details{
Under Brownian motion, the slope of a regression on the returned data frame is expected to be -0.5. Stasis predicts a slope of 
-1, and constant-rate directional evolution a slope of 0.}
\value{
a \code{data.frame} indicating which two tips are being compared, the log10.Rate and the log10.Interval between them}

\references{Gingerich, P. D. 1993 \emph{American Journal of Science} \bold{293A}, 453--478, 
for LRLI plots, though the phylogenetic context is novel. }

\author{Andy Purvis}

\section{Version }{1.3}
\section{Date }{19 May 2006}

\section{Warning }{At present there is no error checking}





## The function is currently defined as
function(phylogeny, phylogeny.matrix, pairs, data.matrix, plot="n"){

#Plots and models LRLI based on phylogenetically independent matched pairs
#'phylogeny' is the phylogeny
#'phylogeny.matrix' is a matrix of PD as constructed by dist.phylo
#'pairs' is an n x 2 matric of tip labels
#'data.matrix' is a matrix of data differences as constructed by as.matrix(dist())
#Returns a data frame (used to be a plot and a model).

#10/5/06: Bug fix to prevent a rate of zero being found (which, when logged, causes error).
#See also all.LRLI

#10/5/06: Bug fix - was not previously passed the phylogeny, which it needed for the tip.labels

#10/5/06: Structure changed to return data frame rather than model.

indices<-matrix(match(pairs[,], phylogeny$tip.label),nrow=length(pairs[,1]),ncol=2)
PD<-rep(NA,length(pairs[,1]))
X.diff<-rep(NA,length(pairs[,1]))

for (i in 1:length(PD))
{
	PD[i]<-phylogeny.matrix[indices[i,1],indices[i,2]]
	X.diff[i]<-data.matrix[indices[i,1],indices[i,2]]
  }

X.diff[X.diff==0]<-min(X.diff[X.diff>0])/10 #Hack
Log10.Rate<-log10(X.diff/PD)
Log10.Interval<-log10(PD)

to.return<-data.frame(pairs[,1],pairs[,2],Log10.Rate,Log10.Interval)
names(to.return)[1:2]<-c("Tip.1","Tip.2")

return(to.return)
  }
}
\keyword{ manip}