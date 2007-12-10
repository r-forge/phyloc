\name{Lacertid.IC}
\alias{Lacertid.IC}
\docType{data}
\title{ Independet contrasts for Bauwens and Diaz-Uriarte (1997) lacertid data set }
\description{
  Independent contrast for Bauwens and Diaz-Uriarte (1997) data set; they are
  based on the data in the file Lacertid.Original using the
  phylogeny in Buwens and Diaz-Uriarte (1997), Tree A  --- see \code{\link{Lacertid.varcov}}}
\format{ The data frame contains eight columns.
  The first seven columns are the (standardized) independent contrasts for the respective variables
  ---see detailed explanation in  \code{\link{Lacertid.Original}}. The
  final column is \describe{
  \item{contr}{the names of the contrasts, as the names of the two
    nodes that form the contrast}}
}




\source{Bauwens, D., and Diaz-Uriarte, R. (1997) Covariation of life-history traits in lacertid lizards: a
  comparative study. \emph{The American Naturalist}, \bold{149}, 91-11}
}


\seealso{ \code{\link{Lacertid.varcov}}, \code{\link{Lacertid.Original}}}
\examples{
# Obtaining correlations through the origin;
# compare with Table 3 in Bauwens and Diaz-Uriarte (1997).

data(Lacertid.IC)
cor.lacert <- matrix(nrow=7,ncol=7) 
for (i in 1:7) for (j in 1:7)
cor.lacert[i,j] <- cor.origin(Lacertid.IC[[i]],Lacertid.IC[[j]])
cor.lacert




\dontrun{
# This data frame can be obtained from the fic data files as:
  
LacertidIC <- cbind(read.table("ifsmi.fic")[,c(3,4)],
                    read.table("ihshw.fic")[,c(3,4)],
                    read.table("iclag.fic")[,c(3,4)],
                    read.table("icfxx.fic")[,3])
stand <- read.table("ifsmi.fic")[,5]
LacertidIC <- LacertidIC/stand
LacertidIC$contr <- read.table("ifsmi.fic")[,1]
names(LacertidIC) <- c("svl","svl.matur", "hatsvl", "hatweight",
                      "clutch.size", "age.mat","cl.freq","contr")
}



}
\keyword{datasets}


