\name{analyse.DI}
\alias{analyse.DI}

\title{ Performs phylogenetic analysis of trait differences vs intervals }
\description{
From a phylogeny and a vector of comparative data, computes, bins and models how trait differences 
scale with time intervals on all independent matched pairs.}
\usage{
analyse.DI(phylogeny, data, log.interval = FALSE, nbins = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ object of class \code{phylo}, assumed to be binary }
  \item{data}{ numeric vector of data for all the terminals in \code{phylogeny}. NOTE: assumed to be
  in same order as \code{phylogeny$tip.label} but this is not checked! A later version will be
  restructured to prevent the obvious scope for terrible mistakes.}
  \item{log.interval}{ Are intervals log-transformed? Default = FALSE. }
  \item{nbins}{ How many bins should comparisons be grouped into? Default is from Sturges' formula. }
}
\details{
  Useful for currently-ongoing investigations into character evolution by Andy Purvis and John Halley
}
\value{
A \code{list} of two linear models, one weighted and one unweighted, of difference on (optionally logged) 
interval. A plot is also produced.}

\author{Andy Purvis}

\section{Version }{1.3}
\section{Date }{19 May 2006}

\section{Warning }{At present there is no error checking}
\section{Warning }{If species values are not in the same order as \code{phylogeny$tip.label}, the function will
return plausible garbage}


## The function is currently defined as
function(phylogeny, data, log.interval=FALSE, nbins=NULL){

#Performs difference vs interval analyses for a phylogeny.

#10/5/06: Produced de novo; links to matched.pairs.DI and binned.DI

phylogeny.mat<-cophenetic(phylogeny)
data.mat<-as.matrix(dist(data))
all.matched.pairs<-find.independent.pairs(phylogeny)

di<-matched.pairs.DI(phylogeny, phylogeny.mat, all.matched.pairs, data.mat)
bdi<-binned.DI(di, log.interval=log.interval, nbins=nbins)

weighted.model<-lm(bdi[,1]~bdi[,2],weights=bdi[,3])
unweighted.model<-lm(bdi[,1]~bdi[,2])
subhead<-paste("Wt slope =",format(as.numeric(coef(weighted.model)[2])),": Unwt slope = ",format(as.numeric(coef(unweighted.model)[2])))
plot(bdi[,1]~bdi[,2],xlab=names(bdi)[2],ylab=names(bdi)[1],cex.lab=1.3,sub=subhead, main=paste("Sample size =",sum(bdi$Sample.Size)))
abline(weighted.model)
abline(unweighted.model,lty=2)

to.return<-list(weighted.model,unweighted.model)

return(to.return)

  }
}
\keyword{ regression }
\keyword{ plot }