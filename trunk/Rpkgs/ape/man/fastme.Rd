\name{FastME}
\alias{fastme}
\alias{fastme.bal}
\alias{fastme.ols}
\title{
  Tree Estimation Based on the Minimum Evolution Algorithm
}
\description{
  The two FastME functions (balanced and OLS) perform the
  Minimum Evolution algorithm of Desper and Gascuel (2002).
}
\usage{
  fastme.bal(X, nni = TRUE)
  fastme.ols(X, nni = TRUE)
}
\arguments{
  \item{X}{a distance matrix; may be an object of class \code{"dist"}.}
  \item{nni}{a boolean value; TRUE to do NNIs (default).}
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Desper, R. and Gascuel, O. (2002) Fast and accurate phylogeny
  reconstruction algorithms based on the minimum-evolution principle.
  \emph{Journal of Computational Biology}, \bold{9(5)}, 687--705.
}
\author{
  original C code by Richard Desper; adapted and ported to R
  by Vincent Lefort \email{vincent.lefort@lirmm.fr}
}
\seealso{
  \code{\link{nj}}, \code{\link{bionj}},
  \code{\link{write.tree}}, \code{\link{read.tree}},
  \code{\link{dist.dna}}, \code{\link{mlphylo}}
}
\examples{
### From Saitou and Nei (1987, Table 1):
x <- c(7, 8, 11, 13, 16, 13, 17, 5, 8, 10, 13,
       10, 14, 5, 7, 10, 7, 11, 8, 11, 8, 12,
       5, 6, 10, 9, 13, 8)
M <- matrix(0, 8, 8)
M[row(M) > col(M)] <- x
M[row(M) < col(M)] <- x
rownames(M) <- colnames(M) <- 1:8
tr <- fastme.bal(M)
plot(tr, "u")
### a less theoretical example
data(woodmouse)
trw <- fastme.bal(dist.dna(woodmouse))
plot(trw)
}
\keyword{models}

