\name{summary.caic}
\alias{summary.caic}
\alias{print.caic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize a caic or macrocaic analysis}
\description{
The summary method simply returns the linear model summary from the caic object. The print method prints some basic information about the analysis followed by the model summary.}
\usage{
\method{summary}{caic}(object, ...)
\method{print}{caic}(x,...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{}
  \item{x}{An object of class 'caic'.}
  \item{\dots}{Arguments to be passed to 'summary.lm'.}
}
\value{
The summary method returns an object of class 'summary.lm'.
}

\author{David Orme}
\seealso{\code{\link{caic}}, \code{link{macrocaic}}}
\examples{

data(SmallTree)
caicMod <- caic(othervar ~ weight, phy=diTree, data=SmallTreeDat, names.col=tip)
print(caicMod)

caicModSummary <- summary(caicMod)
coef(caicModSummary)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
