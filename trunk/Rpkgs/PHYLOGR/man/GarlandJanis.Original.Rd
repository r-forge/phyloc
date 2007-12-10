\name{GarlandJanis.Original}
\alias{GarlandJanis.Original}
\docType{data}
\title{ Garland \& Janis's 1993 data on mammalian running speed and limb length }
\description{
This data set was used by Garland \& Janis in their analysis of
metatarsal/femur ration and running speed in cursorial mammals.  The
data refer to several ecomorphological characteristics for a set of 49
mammals (18 carnivores and 29 ungulates).  
}
\format{
This data frame contains the following columns:
\describe{
\item{Tips}{ the code for each species}
\item{body.mass}{ log 10 of body mass in kilograms}
\item{running.speed}{ log 10 running or sprint speed in km/h}
\item{hind.l.length}{ log 10 hind limb length ---sum of femur, tibia,
  and metatarsal lengths---in cm}
\item{mtf.ratio}{metatarsal/femur ratio}
\item{clade}{ a factor with levels \code{Carnivore} or \code{Herbivore}} 
}
}

  
\source{
Garland, T. Jr., and Janis, C. M. (1993). Does metatarsal/femur ratio
predict maximal running speed in cursorial mammals? \emph{J. Zoology,
  London}, \bold{229}, 133--151.
}

\seealso{ \code{\link{GarlandJanis.IC}}, \code{\link{GarlandJanis.varcov}}}

\examples{
# an example of a GLS fit
data(GarlandJanis.Original)
data(GarlandJanis.varcov)
attach(GarlandJanis.Original)
fit.gls.GJ <- phylog.gls.fit(cbind(body.mass,hind.l.length),
running.speed, GarlandJanis.varcov)
summary(fit.gls.GJ) # summary of the gls model; same as with IC
detach(GarlandJanis.Original)


\dontrun{
# This data set can be prepared from the original pdi files
# (in directory Examples) as:
GarlandJanis.Orig <- read.pdi.data(c("49ms.pdi","49hmt.pdi"),
                   variable.names = c("body.mass", "running.speed",
                                      "hind.l.length","mtf.ratio")) 
Garland.Janis.Orig$clade <- as.factor(c(rep("Carnivore",19),
				      rep("Herbivore",30)))
  }




}
\keyword{datasets}
