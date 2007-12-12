\name{compare.independent.pairs}
\alias{compare.independent.pairs}

\title{ Generate BRUNCH-like contrasts on a fully binary phylogeny}
\description{
Decomposes a fully bifurcating phylogeny into independent pairs of species, which are then 
contrasted}
\usage{
compare.independent.pairs(phylogeny, dataset, variable.list, cherries.only = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ an object of class \code{"phylo"}, assumed to be fully binary. }
  \item{dataset}{ a \code{data.frame} of comparative data, whose first column holds taxon names. }
  \item{variable.list}{ a vector of names of variables to be contrasted }
  \item{cherries.only}{ logical. If TRUE, only cherries are compared; if FALSE, all independent pairs 
  are compared }
}
\details{
In BRUNCH-style contrasts, each species can be used in at most one contrast. No model of evolution 
need be assumed, potentially leading to more robust hypothesis-tests about correlated evolution.}
\value{
A list:
  \item{unscaled.contrasts }{An array of raw differences between compared species, with one column for 
  each element of \code{variable.list}.}
  \item{variances }{The expected variances of the raw differences, under the Brownian motion model of evolution; 
  these may be used for scaling the comparisons if so desired}
  ...
}
\references{ Purvis, A. and Rambaut, A. 1995 \emph{Comput. Appls. Biosci.} \bold{11}, 247--251, for BRUNCH algorithm. }
\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{There is no error-checking}

\seealso{ \code{\link{find.cherries}}, \code{\link{brunch}}}

## The function is currently defined as
function(phylogeny, dataset, variable.list,cherries.only=FALSE){

#Generates BRUNCH-like contrasts in vars from a fully binary phylogeny.

dataset<-strip.spp.not.in.tree(dataset,phylogeny)
dataset<-strip.missing.data(dataset,variable.list)
phylogeny<-strip.phylo(phylogeny,dataset)
nc<-floor(length(dataset[,1])/2)
if (cherries.only==FALSE) pairs<-find.independent.pairs(phylogeny) else pairs<-find.cherries(phylogeny)
names<-list(NULL,variable.list)
contrasts<-array(data=NA,dim=c(length(pairs)/2,length(variable.list)),dimnames=names)
variances<-array(data=NA,dim=length(pairs)/2)

for (i in 1:length(variable.list))
{
	species.data<-dataset[,names(dataset)==variable.list[i]]
	names(species.data)<-dataset[,1]
	pairs.indices<-matrix(match(pairs,names(species.data)),ncol=2)
	contrasts[,i]<-species.data[pairs.indices[,1]]-species.data[pairs.indices[,2]]
  }


phylogeny.mat<-dist.phylo(phylogeny)

for (i in 1:length(contrasts)/2)
{
	variances[i]<-phylogeny.mat[pairs.indices[i,1],pairs.indices[i,2]]
  }


to.return<-list(contrasts,variances)
names(to.return)<-c("unscaled.contrasts","variances")
to.return
  }
}
\keyword{ htest}
\keyword{ manip}