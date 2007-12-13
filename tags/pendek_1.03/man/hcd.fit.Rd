\name{hcd.fit}
\alias{hcd.fit}
\title{Fit an HCD to real data.}
\description{
Takes an integer vector of subclade richness within higher taxa and fits a Hollow Curve distribution to the data. }
\usage{
hcd.fit(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An integer vector of subclade richness values.}
  \item{\dots}{Arguments to be passed to \code{hcd.generate}.}
}
\details{
This function is simply a wrapper to \code{hcd.generate} and takes the number od species and taxa as the sum and length of the provided integer vector. The data is retained in the 'hcd' object produced for use in the plot method.
}
\value{
An object of class 'hcd', as in \code{hcd.generate}, with the addition of two further components:
  \item{real }{The provided integer vector of real taxon richness values}
  \item{data }{The name of the object provided as 'x'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{David Orme}

\seealso{ \code{\link{hcd.generate}} }
\examples{

## 15 Million year old primate taxa from Purvis et al. (1995)
primates <- c(80, 19, 15, 10, 10, 9, 9, 7, 7, 6, 6, 5, 5, 5, 4, 3, 1, 1, 1)
names(primates) <- c("Cercopithecidae","Callitrichidae","Aotus & Callicebus","Pithecia & Cacajao","Galigidae","Hylobates","Lemurinae","Atelini","Cheirogalidae","Aloutta","Indroidea","Pongidae & Homo","Saimiri","Lorisidae","Tarsius","other Cebus","Cebus apella","Lepilemur","Varecia")

hcdPrimate <- hcd.fit(primates)
plot(hcdPrimate)

hcdPrimate <- hcd.fit(primates, minmax=TRUE)
plot(hcdPrimate)


}
\keyword{ distribution }% at least one, from doc/KEYWORDS
