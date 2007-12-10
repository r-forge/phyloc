\name{GarlandJanis.IC}
\alias{GarlandJanis.IC}
\docType{data}
\title{ Independet contrasts for Garland \& Janis data set }
\description{
  Independent contrast for Garland \& Janis (1993) data set; they are
  based on the data in the file GarlandJanis.Original using the
  phylogeny in Garland et al. (1993) --- see GarlandJanis.varcov}
\format{ The data frame contains seven columns.
  The first four columns are the (standardized) independent contrasts for the respective variables
  ---see detailed explanation in
  \code{\link{GarlandJanis.Original}}. The rest of the columns are:
  \describe{
  \item{branch.lengths}{Standard Deviation of Contrast = square root of sum of corrected branch lengths.} 
  \item{names.contr}{the names of the contrasts, as the names of the two
    nodes that form the contrast}
  \item{clade.contr}{ A factor with levels \code{Carnivore},
    \code{Herbivore}, \code{root}, that indicates whether the contrast
    is a contrast within the carnivore lineage, within the ungulates,
    or the root contrast ---between the carnivore and the ungulate clades.}}
}



\source{
Garland, T. Jr., and Janis, C. M. (1993). Does metatarsal/femur ratio
predict maximal running speed in cursorial mammals? \emph{J. zoology,
  London}, \bold{229}, 133--151.

Garland, T. Jr., Dickerman, A. W., Janis, C. M., Jones, J. A. (1993)
Phylogenetic analysis of covariance by computer
simulation. \emph{Systematic Biology}, \bold{42}, 265--292.
}

\seealso{ \code{\link{GarlandJanis.Original}}, \code{\link{GarlandJanis.varcov}}}



\examples{

# Multiple regression with independent contrasts
# excluding the polar bear - grizzly bear contrast
data(GarlandJanis.IC)
lm(running.speed ~ body.mass + hind.l.length - 1,
subset=names.contr!="Tm-Ur", data= GarlandJanis.IC)



\dontrun{
# This data set can be obtained from the original files as:

garland.janis.ic <- cbind(read.table("49ms.fic")[,c(3,4)],
                          read.table("49hmt.fic")[,c(3,4)])

branch.lengths <- read.table("49ms.fic")[,5]
garland.janis.ic <- garland.janis.ic/branch.lengths
names(garland.janis.ic) <- c("body.mass", "running.speed",
                             "hind.l.length", "mtf.ratio")
garland.janis.ic[garland.janis.ic$body.mass<0,] <-
        -1 * garland.janis.ic[ garland.janis.ic$body.mass<0, ]
garland.janis.ic$branch.lengths <- branch.lengths
garland.janis.ic$names.contr <-
                       as.factor(read.table("49ms.fic")[,1])
garland.janis.ic$clade.contr <-
    as.factor( c("root",rep("Carnivore",18), rep("Herbivore",29))) 
}


}
\keyword{datasets}
