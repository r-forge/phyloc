\name{hcd.generate}
\alias{hcd.generate}
\alias{print.hcd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a Hollow Curve Distribution}
\description{
Uses randomization to generate the distribution of subclades among higher taxa under a null model of clade growth.}
\usage{
hcd.generate(species, taxa, reps = 1000, hi = NULL, lo = NULL, minmax = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{ An integer giving the number of species (or more generally any equally ranked subtaxa), to be distributed between higher level taxa. }
  \item{taxa}{An integer giving the number of higher taxa between which to distribute the number of species}
  \item{reps}{The number of randomized curves to distribute}
  \item{hi}{Optional upper limit - the function will record the total number of taxa in the randomizations that greater than or equal to  \code{hi}.}
  \item{lo}{Optional lower limit - the function will record the total number of taxa in the randomizations that less than or equal to  \code{lo}.}
  \item{minmax}{Logical indicating whether or not to keep a record of the minimum and maximum clade size within the ranked taxa. This can be plotted to show an envelope around the average curve.}
}
\details{
Given a number of species drawn from a number of higher level taxa, the number of species within each higher taxon are expected to demonstrate a hollow curve distribution under a null Markov model of clade growth. The HCD is obtained by plotting the species richness of each taxon in descending rank order against their rank. The hollow curve distribution assumes that the higher taxa are of equal age.

This function uses randomization to generate hollow curve distributions showing the expectations under this null model. The function simulates a broken stick procedure and returns the average species richness at each rank in the HCD.
}
\value{
An object of class 'hcd', containing the following components:
\item{species}{The integer provided as \code{species}.}
\item{taxa}{The integer provided as \code{taxa}.}
\item{reps}{The number of replicates used to generate the HCD.}
\item{means}{A vector giving the average species richness in each taxon, ranked by species richness.}
\item{min}{If \code{minmax} is TRUE then the function returns a vector giving the least species rich taxon at each rank.}
\item{max}{If \code{minmax} is TRUE then the function returns a vector giving the most species rich taxon at each rank.}
\item{hi}{If a value is provided for \code{hi} then that is recorded here.}
\item{num.hi}{If a value is provided for \code{hi} then \code{num.hi} gives the total number of taxa across all simulations that are greater than or equal to \code{hi}.}
\item{lo}{If a value is provided for \code{lo} then that is recorded here.}
\item{num.lo}{If a value is provided for \code{lo} then \code{num.lo} gives the total number of taxa across all simulations that are less than or equal to \code{lo}.}
}
\references{  }
\author{David Orme}

\seealso{\code{\link{hcd.fit}}, \code{\link{export.hcd}} }
\examples{

hcdExample <- hcd.generate(species=50, taxa=20)
print(hcdExample)
plot(hcdExample)

hcdExample <- hcd.generate(species=50, taxa=20, lo=2, hi=5, minmax=TRUE)
print(hcdExample)
plot(hcdExample)


}
\keyword{distribution}
