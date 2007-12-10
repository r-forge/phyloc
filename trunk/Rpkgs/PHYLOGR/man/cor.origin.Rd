\name{cor.origin}
\alias{cor.origin}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Correlation Through the Origin}
\description{
 Return the correlation through the origin of two vectors. Generally
 used for indepdendent contrasts
}
\usage{
cor.origin(x, y)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ A vector }
  \item{y}{ A vector (of same size as x) }
}
\value{
  The correlation of x and y, from a model without intercept (i.e.,
  forcing the line through the origin).
}
\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
}
\author{ R. Diaz-Uriarte and T. Garland, Jr.}
\note{ This is a very simple function, provided for convenience.  You
  can obtain the p-value, if you wish, with the usual formula for the
  t-statistic:
  \code{2*(1 - pt(sqrt(df) * rho / sqrt(1 - rho^2), df))}
  where rho is the correlation through the origin and df are the appropriate degrees of
freedom ---generally N-1.
}


\examples{
 x <- rnorm(100)
 y <- rnorm(100)
 rho <- cor.origin(x,y)
 rho # the correlation
 2 * (1 - pt(sqrt(99) * rho / sqrt(1 - rho^2), 99))  # the p-value
 
}
\keyword{ models }%-- one or more ...
