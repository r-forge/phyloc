\name{export.hcd}
\alias{export.hcd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Export a hollow curve distribution.}
\description{
Saves an object of class 'hcd' to a text file.}
\usage{
export.hcd(x, outfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class 'hcd'.}
  \item{outfile}{The name of the text file to save to.}
}
\details{
The file consists of a number of commented out header lines describing how the data were generated, followed by a tab delimited table of the means, along with real richness and minmax values if present in the 'hcd' object.
}
\author{David Orme}
\seealso{\code{\link{hcd.fit}}, \code{\link{hcd.generate}}}
\examples{

primates <- c(80, 19, 15, 10, 10, 9, 9, 7, 7, 6, 6, 5, 5, 5, 4, 3, 1, 1, 1)
names(primates) <- c("Cercopithecidae","Callitrichidae","Aotus & Callicebus","Pithecia & Cacajao","Galigidae","Hylobates","Lemurinae","Atelini","Cheirogalidae","Aloutta","Indroidea","Pongidae & Homo","Saimiri","Lorisidae","Tarsius","other Cebus","Cebus apella","Lepilemur","Varecia")

hcdPrimate <- hcd.fit(primates, minmax=TRUE)
export.hcd(hcdPrimate, "primateHCD.txt")


}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ distribution }% __ONLY ONE__ keyword per line
