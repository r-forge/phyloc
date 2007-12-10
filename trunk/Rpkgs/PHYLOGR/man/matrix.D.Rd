\name{matrix.D}
\alias{matrix.D}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Produce the Matrix D of Garland and Ives (2000) for GLS}
\description{
  Produces the matrix D in Garland and Ives, 2000, p. 361, for further use
  in GLS procedures. You will rarely need to call this function
  directly. It is called by phylog.gls.fit.
}
\usage{
matrix.D(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ a phylogenetic variance-covariance matrix, such as a matrix returned
   from a call to read.phylog.matrix. }
 
}
\value{A variance-covariance matrix.}
\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
  
  Garland, T. Jr. and Ives, A. R. (2000) Using the past to predict the
  present: confidence intervals for regression equations in phylogenetic
  comparative methods. \emph{The American Naturalist}, 155, 346-364.
}

\note{ The file read is a dsc file generated from PDDIST under option 5
  with the additional options : - in matrix form; - with header; - with
  scaled values. }
\author{Ramon Diaz-Uriarte and Theodore Garland, Jr.}


\seealso{\code{\link{read.phylog.matrix}}, \code{\link{phylog.gls.fit}}}

\examples{
#perform GLS using lm function after transforming the variables
data(Lacertid.varcov)
mD <- matrix.D(Lacertid.varcov)
data(Lacertid.Original)
# obtain transformed variables
Z <- mD \%*\% Lacertid.Original$clutch.size
U <- mD \%*\% cbind(rep(1,18),Lacertid.Original$svl) # intercept included
lm1 <- lm(Z ~ U - 1) # eliminate intercept, since already included in U matrix
summary(lm1)
}

\keyword{array}%-- one or more ...
