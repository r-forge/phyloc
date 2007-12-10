\name{GarlandJanis.varcov}
\alias{GarlandJanis.varcov}
\title{ Phylogenetic variance-covariance matrix for Garland \& Janis (1993). }
\description{
  Phylogenetic variance-covariance matrix for the species in Garland \&
  Janis (1993). Note that the phylogeny is not exactly the same as in
  Garland \& Janis (1993), but actually corresponds to the more recent
  phylogeny in Garland et al. (1993).
}
\format{ A matrix with the phylogenetic distance between species; every
  entry dij is the sum of branch segment lengths that species i and j share in common.}
\source{
Garland, T. Jr., and Janis, C. M. (1993). Does metatarsal/femur ratio
predict maximal running speed in cursorial mammals? \emph{J. Zoology,
  London}, \bold{229}, 133--151.

Garland, T. Jr., Dickerman, A. W., Janis, C. M., Jones, J. A. (1993)
Phylogenetic analysis of covariance by computer
simulation. \emph{Systematic Biology}, \bold{42}, 265--292.
}

\seealso{ \code{\link{GarlandJanis.Original}}, \code{\link{GarlandJanis.IC}}}

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
# This data set can be obtained from the original dsc file as:
  
    GarlandJanis.varcov <- read.phylog.matrix("49ms.dsc")


}



}


\keyword{datasets}
