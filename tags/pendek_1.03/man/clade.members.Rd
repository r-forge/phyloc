\name{clade.members}
\alias{clade.members}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Obtain a vector of the tips subtending from a node in a phylogeny}
\description{
   Returns a vector of the tips that subtend from a supplied node in a phylogeny.
}
\usage{
clade.members(x, phyl, tip.labels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A integer in the range minus number of internal of nodes to positive number of tips. }
  \item{phyl}{An object of class 'phylo'. }
  \item{tip.labels}{A logical flag indicating whether to return the node numbers of the tips or their tip labels.}
}
\details{

}
\value{
Either a numeric vector of external node (i.e. tip) numbers or a character vector of tip labels.
}
\author{David Orme}


\seealso{\code{\link{all.clades}}, \code{\link{clade.stats}} }
\examples{
data(bird.orders)
clade.members(-3, bird.orders, tip.labels=TRUE)
}
\keyword{ manip }% at least one, from doc/KEYWORDS
\keyword{ utilities}% __ONLY ONE__ keyword per line
