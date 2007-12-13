\name{test.all.linear}
\alias{test.all.linear}
\title{ Test linearity of multiple single-predictor relationships}
\description{
Tests for significance of cubic and quadratic terms in comparative relationships between each of 
several predictor variables in turn and a single response variable
}
\usage{
test.all.linear(phylogeny, dataset, response, predictors, check.robust = FALSE, cutoff = 3, use.robust = FALSE, kappa = NULL, filter = NULL)
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
  \item{filter}{Vector of names holding the names of further columns holding dummy variables to permit 
  systematic exclusion of cases: taxa for which these columns contain missing values will be excluded from 
  the analysis. }
}
\details{
  Basically just a wrapper for \code{test.linear}
}
\value{
Nothing is returned at present; all output goes to the screen.}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking}

\seealso{ \code{\link{test.linear}},  \code{\link{pic.lm}}}


## The function is currently defined as
function(phylogeny, dataset, response, predictors, check.robust=FALSE, cutoff=3, use.robust=FALSE, kappa=NULL, filter=NULL){

for(i in 1:length(predictors))
{
	print(paste("Computing model for",predictors[i]))
	test.linear(phylogeny, dataset, c(response,predictors[i]),check.robust=check.robust, cutoff=cutoff, use.robust=use.robust,kappa=kappa,filter=filter)
  }

  }
}
\keyword{ regression }
\keyword{ nonlinear }
