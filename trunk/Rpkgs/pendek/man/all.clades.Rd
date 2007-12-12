\name{all.clades}
\alias{all.clades}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Obtain a list of members for internal nodes}
\description{
Returns a list of vectors showing the tips that descend from each node in a phylogeny.}
\usage{
all.clades(phyl, tips = FALSE, tip.labels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phyl}{ An object of class 'phylo'. }
  \item{tips}{ A logical indicating whether to include external node membership in the list.}
  \item{tip.labels}{A logical indicating whether to return numeric tip codes or character tip.labels.}
}
\details{
The function runs \code{\link{clade.members}} over each node in the phylogeny, possibly including the external nodes as indicated by the tips argument, and returns a list of vectors showing the members of the clade defined by each node.}
\value{
A list of either integers or characters where each list component shows the tips  descending from each node. The components are labeled by node number.
}
\author{David Orme}
\seealso{\code{\link{clade.members}}}
\examples{
data(bird.orders)
all.clades(bird.orders)
}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{utilities}% __ONLY ONE__ keyword per line
