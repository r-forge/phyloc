\name{optimal.pic}
\alias{optimal.pic}
\alias{optimise.kappa}

\title{Phylongenetically independent contrasts with kappa optimisation}
\description{
Takes a variable and a phylogeny, finds the optimal value of kappa (\code{optimise.kappa}), transforms the 
branch lengths accordingly, and computes and returns phylogenetically independent contrasts}
\usage{
optimal.pic(phylogeny, dataset, variable)
optimise.kappa(phylogeny, dataset, variable)
}

\arguments{
  \item{phylogeny}{ An object of class \code{"phylo"} which must not contain polytomies }
  \item{dataset}{ A data frame whose first column holds taxon names }
  \item{variable}{ Name of column holding the variable in which contrasts are required }
}
\details{
\code{optimal.pic} is used by \code{pic.lm} and associated functions, and calls 
\code{optimise.kappa}. Kappa is a branch length scaling parameter, optimised by minimising 
the correlation between scaled contrasts and the scaling factor (square root of total branch length).}
\value{
\code{optimal.pic} returns a numeric vector of phylogenetically independent contrasts;
\code{optimise.kappa} returns a single number.}
\references{Garland, T. Jr., Harvey, P. H., and Ives, A.R. 1992 \emph{Systematic 
Biology} \bold{41}, 18--32}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking}

\seealso{ \code{\link{pic.lm}}, \code{\link{optimise.kappa}}}

## The function is currently defined as
function(phylogeny, dataset, variable){

#Written by Andy Purvis
#Returns set of contrasts computed with the optimally transformed branch lengths
#variable is the name of the variable of interest

#requires ape

kappa<-optimise.kappa(phylogeny, dataset, variable)
phylogeny$edge.length<-phylogeny$edge.length^kappa
species.data<-dataset[,names(dataset)==variable]
names(species.data)<-dataset[,1]
dx<-pic(species.data,phylogeny)
dx
  }
}

\keyword{regression}