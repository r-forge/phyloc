\name{Lacertid.varcov}
\alias{Lacertid.varcov}
\docType{data}
\title{ Variance-covariance matrix for lacertids from Bauwens and
  Diaz-uriarte (1997) }
\description{ Phylogenetic variance-covariance matrix for 18 species of
  lacertids. It is based on Tree A of Bauwens and Diaz-Uriarte (1997).}

\format{
The \code{Lacertid.varcov} data frame has 18 rows and 18 columns,
corresponding to each one of the 18 lacertidspecies; the matrix is the
phylogenetic variance-covariance matrix between all 18 species; thus,
each entry dij is the the sum of branch segment lengths that species i and j share in common.

}

  
\source{Bauwens, D., and Diaz-Uriarte, R. (1997) Covariation of life-history traits in lacertid lizards: a
  comparative study. \emph{The American Naturalist}, \bold{149}, 91-11}
}


\seealso{ \code{\link{SimulExample}}, \code{\link{Lacertid.Original}}}
\examples{
# a GLS fit
data(Lacertid.varcov)
data(Lacertid.Original)
ex.gls.phylog <-
phylog.gls.fit(Lacertid.Original$svl,Lacertid.Original$clutch.size,Lacertid.varcov)
ex.gls.phylog


\dontrun{
# This data can be obtained from the original dsc file as:
Lacertid.varcov <- read.phylog.matrix("ifsmi.dsc")
  }



}
\keyword{datasets}
