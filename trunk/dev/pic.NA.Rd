\name{pic.NA}
\alias{pic.NA}
\title{Phylogenetically Independent Contrasts with missing data}
\usage{
pic.NA(x, phy, ...)
}
\arguments{
  \item{x}{a numeric data.frame of data, may contain missing values.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{scaled}{optional arguments to be passed to \code{"pic"}}
}
\description{
  Compute phylogenetically independent contrasts for multiple traits with missing values.
}
\details{
  This function is a wrapper for the \code{\link{pic}} function. It calculates phylogenetically independent contrasts for a data.frame. Each column of the data frame contains trait data, which may include missing values. For each trait, the \code{"phylo"} object is pruned to include only the taxa with non-missing data for that trait and contrasts are calculated.
  
  If \code{x} has row.names, its values are matched to the tip labels of
  \code{phy}, otherwise its values are taken to be in the same order
  than the tip labels of \code{phy}.

  The user must be careful here since the function requires that both
  series of names perfectly match, so this operation may fail if there
  is a typing or syntax error. If both series of names do not match, the
  values in the \code{x} are taken to be in the same order than the tip
  labels of \code{phy}, and a warning message is issued.
}
\value{
  A \code{data.frame} containing contrasts for each node in the \code{phylo}. Columns in the data frame correspond to the columns of \code{x}. NA values indicate contrasts that were not calculated since the corresponding node was pruned due to missing trait values at the tips.
}
\references{
  Felsenstein, J. (1985) Phylogenies and the comparative method.
  \emph{American Naturalist}, \bold{125}, 1--15.
}
\author{Steven Kembel \email{skembel@berkeley.edu}}
\seealso{
  \code{\link{pic}}
}
\examples{
tr <- rcoal(20)
dat <- replicate(10,evolve.brownian(tr))
dat[1:5,1:5] <- NA
pic.missing(dat,tr)
}
\keyword{regression}
