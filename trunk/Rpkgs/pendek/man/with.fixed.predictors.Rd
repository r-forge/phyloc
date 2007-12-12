\name{with.fixed.predictors}
\alias{with.fixed.predictors}

\title{ Fit multiple pic.lm models in which the response and some predictors are common to all }
\description{
This function is a tool for heuristic search of model space for the best minimum adequate model in 
phylogenetic comparative analyses. After \code{pic.mod.simp} has performed its greedy search, it is 
sensible to run this function with the predictors in the model returned by \code{pic.mod.simp} as 
\code{fixed} and those not in the model as \code{predictors}, to see whether any of them joins the 
model. }
\usage{
with.fixed.predictors(phylogeny, dataset, response, fixed, predictors, check.robust = FALSE, cutoff = 3, use.robust = FALSE, kappa = NULL, filter = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{"phylo"}, which can contain polytomies}
  \item{dataset}{A data frame whose first column holds taxon names }
  \item{response}{Name of the column holding the response variable in the model to be fitted }
  \item{fixed}{Names of the columns holding predictors to be included in all models. Typically 
  these will be the variables that were predictors in the model returned by \code{pic.mod/simp}.}
  \item{predictors}{Names of the columns holding other other predictors to be included one at a 
  time in the model.}
  \item{check.robust}{ If \code{TRUE}, re-fits the model excluding points with studentised 
  residuals whose absolute value exceeds \code{cutoff} }
  \item{cutoff}{Number indicating how large a studentised residual must be before a point is excluded }
  \item{use.robust}{If \code{TRUE}, the returned model is that found with \code{check.robust}; by 
  default, the model returned is based on all points}
  \item{kappa}{Optional branch length transformation. If \code{kappa} is specified, branch lengths are 
  raised to the power of \code{kappa}. If \code{kappa} is not specified, it is optimised for each 
  variable in turn. }
  \item{filter}{Vector of names holding the names of further columns holding dummy variables to permit 
  systematic exclusion of cases: taxa for which these columns contain missing values will be excluded from 
  the analysis. }
}

}
\value{
a \code{list} of models, each of class \code{"lm"}, one for each of \code{predictors}.}

\references{Purvis, A., Gittleman, J. L., Cowlishaw, G. and Mace, G. M. 2000 \emph{Proc. R. Soc. Lond. B} 
\bold{267}, 1947-1952, for the approach}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking}

\seealso{ \code{\link{picmodsimp}}}


## The function is currently defined as
function(phylogeny, dataset, response, fixed, predictors, check.robust=FALSE, cutoff=3, use.robust=FALSE, kappa=NULL, filter=NULL){

models<-as.list(predictors)

for(i in 1:length(predictors))
{
	print(paste("Computing model for",predictors[i]))
	models[[i]]<-summary.pic.lm(pic.lm(phylogeny, dataset, response, c(predictors[i],fixed),check.robust=check.robust, cutoff=cutoff, use.robust=use.robust, kappa=kappa))
  }

models
  }
}
\keyword{regression}