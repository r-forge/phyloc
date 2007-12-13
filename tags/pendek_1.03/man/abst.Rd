\name{abst}
\alias{abst}

\title{Internal helper function for optimising kappa}
\description{
Finds absolute value of t for regression of absolute scaled contrasts on scaling factor.}
\usage{
abst(kappa, phylogeny, dataset, variable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{ numeric. Branch lengths are raised to the power \code{kappa}. }
  \item{phylogeny}{ An object of class \code{"phylo"} which must not contain polytomies }
  \item{dataset}{ A data frame whose first column holds taxon names }
  \item{variable}{ Name of column holding the variable in which contrasts are required }
}

}
\value{
  numeric. A positivised t-statistic.
}
\references{Garland, T. Jr., Harvey, P. H., and Ives, A.R. 1992 \emph{Systematic 
Biology} \bold{41}, 18--32}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking}

\seealso{ \code{\link{optimise.kappa}}}

## The function is currently defined as
function(kappa, phylogeny, dataset, variable){

#Written by Andy Purvis
#Finds absolute value of t for regression of absolute scaled contrasts on scaling factor (SD)

species.data<-dataset[,names(dataset)==variable]
names(species.data)<-dataset[,1]
phylogeny$edge.length<-phylogeny$edge.length^kappa
dx<-pic(species.data,phylogeny,var.contrasts=TRUE)
SD<-sqrt(dx[,2])
model1<-lm(abs(dx[,1])~SD)
T<-summary(model1)$coefficient[2,3]
abst<-abs(T)

abst
  }
}
\keyword{ manip}
\keyword{ regression}