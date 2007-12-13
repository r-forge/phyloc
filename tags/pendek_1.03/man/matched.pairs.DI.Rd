\name{matched.pairs.DI}
\alias{matched.pairs.DI}

\title{Differences and intervals on phylogenetically matched pairs}
\description{
Computes trait differences and time intervals for each of a number of phylogenetically independent 
matched pairs}
\usage{
matched.pairs.DI(phylogeny, phylogeny.matrix, pairs, data.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ An object of class \code{phylo} }
  \item{phylogeny.matrix}{ A matrix of phylogenetic distances between pairs of terminals, as 
  constructed by \code{cophenetic}.  }
  \item{pairs}{ n x 2 matrix of tip labels, as produced by \code{find.independent.pairs} }
  \item{data.matrix}{ matrix of data differences, as constructed by \code{as.matrix(dist())}.  }
}
\details{
  Useful for currently-ongoing investigations into character evolution by Andy Purvis and John Halley
}
\value{
a \code{data.frame} indicating which two tips are being compared, the log10.Rate and the log10.Interval between them}

\author{Andy Purvis}

\section{Version }{1.3}
\section{Date }{19 May 2006}

\section{Warning }{At present there is no error checking}



\seealso{ \code{\link{binned.DI}}, \code{\link{analyse.DI}}}


## The function is currently defined as
function(phylogeny, phylogeny.matrix, pairs, data.matrix){

#Computes differences in a trait and intervals for phylogenetically independent matched pairs
#'phylogeny' is the phylogeny
#'phylogeny.matrix' is a matrix of PD as constructed by cophenetic
#'pairs' is an n x 2 matrix of tip labels
#'data.matrix' is a matrix of data differences as constructed by as.matrix(dist())
#Returns a data frame with four columns, see code for details.


#10/5/06: Created from matched.pairs.LRLI, after bug fixes to that function

indices<-matrix(match(pairs[,], phylogeny$tip.label),nrow=length(pairs[,1]),ncol=2)
PD<-rep(NA,length(pairs[,1]))
X.diff<-rep(NA,length(pairs[,1]))

for (i in 1:length(PD))
{
	PD[i]<-phylogeny.matrix[indices[i,1],indices[i,2]]
	X.diff[i]<-data.matrix[indices[i,1],indices[i,2]]
  }

X.diff[X.diff==0]<-min(X.diff[X.diff>0])/10 #Hack

to.return<-data.frame(pairs[,1],pairs[,2],X.diff,PD)
names(to.return)<-c("Tip.1","Tip.2","Difference","Interval")


return(to.return)
  }
}
\keyword{ manip}