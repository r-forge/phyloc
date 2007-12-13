\name{all.one.predictors}
\alias{all.one.predictors}

\title{ Automatically perform pic.lm for each of several predictors and a fixed response variable }
\description{
Automatically use \code{pic.lm} to regress a given response variable on each of multiple 
predictor variables in turn using phylogenetically independent contrasts}
\usage{
all.one.predictors(phylogeny, dataset, response, predictors, check.robust = FALSE, cutoff = 3, use.robust = FALSE, kappa = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{"phylo"}, which can contain polytomies}
  \item{dataset}{A data frame whose first column holds taxon names }
  \item{response}{Name of the column holding the response variable in the model to be fitted }
  \item{predictors}{ Vector of names of columns holding explanatory variables for the model}
  \item{check.robust}{ If \code{TRUE}, re-fits the model excluding points with studentised 
  residuals whose absolute value exceeds \code{cutoff} }
  \item{cutoff}{Number indicating how large a studentised residual must be before a point is excluded }
  \item{use.robust}{If \code{TRUE}, the returned model is that found with \code{check.robust}; by 
  default, the model returned is based on all points}
  \item{kappa}{Optional branch length transformation. If \code{kappa} is specified, branch lengths are 
  raised to the power of \code{kappa}. If \code{kappa} is not specified, it is optimised for each 
  variable in turn. }
}

}
\value{
A list of regression models, each with class \code{"lm"}, each produced by \code{pic.lm}}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking on the phylogeny passed to the function}

\seealso{ \code{\link{pic.lm}} for more details }


## The function is currently defined as
function(phylogeny, dataset, response, predictors,check.robust=FALSE, cutoff=3, use.robust=FALSE, kappa=NULL){

models<-as.list(predictors)

for(i in 1:length(predictors))
{
	print(paste("Computing model for",predictors[i]))
	models[[i]]<-summary.pic.lm(pic.lm(phylogeny, dataset, response, predictors[i],check.robust=check.robust, cutoff=cutoff, use.robust=use.robust, kappa=kappa))
  }

models
  }
}
\keyword{regression}