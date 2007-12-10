\name{phylog.gls.fit}
\alias{phylog.gls.fit}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Phylogenetically-Based GLS Model Fitting }
\description{
  Fits a GLS linear model, as in Garland and Ives (2000), using a
  phylogenetic variance-covariance function}
\usage{
phylog.gls.fit(x, y, cov.matrix, intercept = TRUE, exclude.tips = NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{The predictor or ''X'' variables; they must be numeric
   variables. If you are using a factor, you must recode it numerically,
 with the appropriate type of contrasts ---see function contrasts and
 Venables & Ripley (1999) ch. 6}
 \item{y}{ The response}
 \item{cov.matrix}{ The phylogenetic variance-covariance matrix, which
   can be obtained from read.phylog.matrix.}
 \item{intercept}{ Include an intercept in the model? Defaults to TRUE }
 \item{exclude.tips}{ The tips to be excluded from the
   analyses. Defaults to NULL}
}

\value{a fitted linear model}
\references{

  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
  
  Garland, T. Jr. and Ives, A. R. (2000) Using the past to predict the
  present: confidence intervals for regression equations in phylogenetic
  comparative methods. \emph{The American Naturalist}, 155, 346-364.

  Venables, W. N. and Ripley, B. D. (1999) \emph{Modern applied
    statistics with S-Plus, 3rd ed.} Springer-Verlag.

}
\author{Ramon Diaz-Uriarte and Theodore Garland, Jr.}

\section{WARNING}{This is one possible implementation of GLS that uses
  the transformation of the Y and X as explained in Garland and Ives (2000). Ideally, we
  would directly call gls from the NLME package, passing it the var-cov
  matrix, but there are some printing problems of the fitted object in the R implementation when we use a
  fixed correation structure. The advantage
  of using gls from NLME is that the function is called using the
  typical syntax for linear models, and we do not need to worry about making
  categorical factors into numerical variables. Once the problem with
  NLME is solved I'll add functions to incorporate GLS into the analysis
  of data sets.
  
  In the meantime, when using this function, you should be aware
  that:
  
  1) the overall F-test reported is wrong (it is like comparing to a
  model without an intercept);
  
  2) you can apply the usual plot(fitted.model) to see diagnostic plots,
  or other diagnostic functions such as lm.influence,
  influence.measures, etc.  But most of these will be wrong and meaningless.}




\seealso{\code{\link{read.phylog.matrix}}, \code{\link{matrix.D}}}

\examples{
data(Lacertid.varcov)
data(Lacertid.Original)
ex.gls.phylog <-
phylog.gls.fit(Lacertid.Original$svl,Lacertid.Original$clutch.size,Lacertid.varcov)
ex.gls.phylog
}
\keyword{models }%-- one or more ...
