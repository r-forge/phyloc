\name{SimulExample}
\alias{SimulExample}
\docType{data}
\title{ A simulated data set }
\description{
  A simulated data set; the phylogeny is based in Bauwens and Diaz-Uriarte
  (1997), such as is included in the file ifsm.pdi (in the Examples
  directory). But the data are all completely fictitious and have
  nothing to do with lacertids (or, for that matter, with any other creatures).}
}
\format{
This data frame contains the following columns:
\describe{
\item{sim.counter}{the simulation counter}
\item{Tips}{ the name of tips; it matches those for the lacertid
  examples but, again, is unrelated to those}
\item{y}{ one numeric variable}
\item{x1}{another numeric variable}
\item{x2}{ditto}
\item{x3}{ditto}
\item{x4}{ditto}
\item{x5}{guess what? same thing}
\item{x6}{ again}
\item{x7}{ once more}
\item{diet}{a factor with fictitious levels
\code{Carnivore} 
\code{Herbivore} 
\code{Ommnivore} 
}
}
}
\source{Bauwens, D., and Diaz-Uriarte, R. (1997) Covariation of life-history traits in lacertid lizards: a
  comparative study. \emph{The American Naturalist}, \bold{149}, 91-11}
}

\examples{
# a canonical correlation example
data(SimulExample)
ex1.cancor <- cancor.phylog(SimulExample[,c(1,2,3,4,5)],SimulExample[,c(1,2,6,7,8)])
ex1.cancor
summary(ex1.cancor)
plot(ex1.cancor)
}
\keyword{datasets}
