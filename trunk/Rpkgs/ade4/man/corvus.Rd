\name{corvus}
\alias{corvus}
\docType{data}
\title{Corvus morphology}
\description{
This data set gives a morphological description of 28 species of the genus Corvus split in two habitat types and phylogeographic stocks.
}
\usage{data(corvus)}
\format{
\code{corvus} is data frame with 28 observations (the species) and 4 variables : 
   \item{wing}{: wing length (cm)}
   \item{bill}{: bill length (cm)}
   \item{habitat}{: habitat with two levels \code{clos} and \code{open}}
   \item{phylog}{: phylogeographic stock with three levels \code{amer}(America), \code{orien}(Oriental-Australian),
    \code{pale}(Paleoarctic-African)}
}
\references{
Laiolo, P. and Rolando, A. (2003) The evolution of vocalisations in the genus Corvus: effects of phylogeny, morphology and habitat. 
\emph{Evolutionary Ecology}, \bold{17}, 111--123.
}
\examples{
data(corvus)
plot(corvus[,1:2])
s.class(corvus[,1:2], corvus[,4]:corvus[,3], add.p = TRUE)
}
\keyword{datasets}
