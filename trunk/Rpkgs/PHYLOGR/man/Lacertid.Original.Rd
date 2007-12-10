\name{Lacertid.Original}
\alias{Lacertid.Original}
\docType{data}
\title{ Bauwens and Diaz-Uriarte (1997) lacertid data }
\description{
This is part of the data set used by Bauwens and Diaz-Uriarte (1997) in
their analysis of lacertid life histories.  The data include several
life history traits of 18 lacertid species.  
}
\format{
This data frame contains the following columns:
\describe{
\item{Tips}{ the code for each species}
\item{svl}{ log 10 of mean adult female Snout-to-Vent length in mm}
\item{svl.matur}{ log 10 of SVL when sexual maturity (females) is reached}
\item{hatsvl}{ log 10 of hatchling svl in mm}
\item{hatweight}{log10 of hatchling mass in grams}
\item{clutch.size}{ log10 of clutch size}
\item{age.mat}{log10 of age at maturity in months}
\item{cl.freq}{ log10 of clutch frequency ---number of clutches per year}
}}

  


\source{Bauwens, D., and Diaz-Uriarte, R. (1997) Covariation of life-history traits in lacertid lizards: a
  comparative study. \emph{The American Naturalist}, \bold{149}, 91-11}
}


\seealso{ \code{\link{Lacertid.varcov}}, \code{\link{Lacertid.IC}}}
\examples{
# a GLS fit
data(Lacertid.varcov)
data(Lacertid.Original)
ex.gls.phylog <-
phylog.gls.fit(Lacertid.Original$svl,Lacertid.Original$clutch.size,Lacertid.varcov)
ex.gls.phylog



\dontrun{
  # This data set can also be obtained from the pdi files
  # (see example in GarlandJanis.Original), or as:

LacertidSim <- read.sim.data(c("ifsmi.sim","ihshw.sim","iclag.sim","icfxx.sim"),
                  pdi.files=c("ifsmi.pdi","ihshw.pdi","iclag.pdi", "icfxx.pdi"),
		  variable.names = c("svl","svl.matur","hatsvl","hatweight",
		                    "clutch.size", "age.mat","cl.freq", "xx"))

LacertidSim <- LacertidSim[,-10]
LacertidOriginal <- LacertidSim[LacertidSim$sim.counter==0,-1]
  }



}
\keyword{datasets}


