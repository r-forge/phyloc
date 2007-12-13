\name{test.linear}
\alias{test.linear}

\title{Test linearity of comparative relationship between two variables}
\description{
Tests whether a comparative relationship between a single predictor variable and a 
single response variable is linear, or has significant quadratic or cubic terms.
}
\usage{
test.linear(phylogeny, dataset, variable.list, check.robust = FALSE, cutoff = 3, use.robust = FALSE, kappa = NULL, filter = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{"phylo"}, which can contain polytomies}
  \item{dataset}{A data frame whose first column holds taxon names }
  \item{variable.list}{A vector holding the names of two variables; the first is taken to be 
  the response variable, and the second the predictor}
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
The predictor is squared and cubed, and contrasts computed in these as well as in the original 
predictor; the three variables are then passed as predictors to \code{pic.lm}}
\value{
At present, nothing is returned; all output is printed to the screen}

\author{Andy Purvis}

\section{Version }{1.3}
\section{Date }{19 May 2006}

\section{Warning }{At present there is no error checking at all}

\seealso{ \code{\link{pic.lm}}; \code{\link{test.all.linear}} for testing linearity of each of 
several predictors in turn}

\keyword{regression}
\keyword{nonlinear}