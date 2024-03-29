\name{name.check}
\alias{name.check}
\alias{treedata}
\title{ Compares taxa in data and tree }
\description{
This function is a general tool for checking for concordance between a data file and a phylogenetic
tree.  For the data, names can be specified as the names of objects in the vector, rownames of the data 
array or as 'data.names'. The name.check function finds and lists all taxa present in data set but not 
in the tree, and vice-versa.  The treedata function returns a list containing both the tree and the data 
after pruning out any species that are not found in both. }
}
\usage{
name.check(phy, data, data.names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ an object of class "phylo" }
  \item{data}{ data for tips of the tree }
  \item{data.names} { names of the tips in the order of the data; if this is not given, names 
  					will be taken from the names or rownames of the object data}
}
\details{
}
\value{
  \item{Tree.not.data}{Taxa in tree but not data}
  \item{Data.not.tree}{Taxa in data but not tree}
  ...
}
\references{  }
\author{ Luke J. Harmon }
\note{ 

}
\seealso{  }
\examples{

data(geospiza)
attach(geospiza)

name.check(geospiza.tree, geospiza.data)

}

\keyword{ arith }
