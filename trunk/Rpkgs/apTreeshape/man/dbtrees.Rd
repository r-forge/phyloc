\name{dbtrees}
\alias{dbtrees}

\title{ Download trees from internet databases.}

\description{
\code{dbtrees} connects to the internet phylogenetic tree database TreeBASE or Pandit, downloads trees with specified ID numbers, and converts them into either the \code{"treeshape"} or the \code{"phylo"} format. Can solve polytomies. 
}

\usage{
dbtrees(db, tree, class="phylo", type="s", quiet=FALSE, 
      model=NULL, p=0.3)
}

\arguments{
\item{db}{ A character string equals to the name of the database to connect: \code{"pandit"} or \code{"treebase"}.}
\item{tree}{ A vector or list of ID corresponding to the trees to be downloaded.}
\item{class}{The class of the downloaded objects: must be one of \code{"phylo"} (default) or \code{"treebase"}.}
\item{type}{ The type of the tree: an object of class \code{"character"} equals to either \code{"s"} (default) or \code{"f"} (seed or full). This option is only used for the Pandit database.}
\item{quiet}{ A \code{logical} value. If \code{TRUE}, nothing will be printed on screen. If \code{FALSE} (default), informations about the status of the download are displayed.}
\item{model}{ Argument to be passed to \code{as.treshape.phylo}. One of \code{NULL} (default), \code{"pda"}, \code{"biased"}, \code{"aldous"} or \code{"yule"}}
\item{p}{ Argument to be passed to \code{as.treshape.phylo} (only for the biased speciation model).}  
}

\details{
The aim of this function is to provide a simple method to retrieve trees and to put them into the \code{"treeshape"} or \code{"phylo"} format, in order to be able to apply the methods of \code{ape} or \code{apTreeshape}. \cr
This function connects to internet database and download trees using their internal reference number. There is currently no option to retrieve trees from names. See references for more details about the Pandit and TreeBASE databases and the description of their data. \cr
This function can solve trees with polytomies using simulation methods (see \code{as.treeshape} for more details about the polytomy resolution).
}

\value{
An object of class \code{"treeshape"} or class \code{"phylo"} if only one tree is requested, a list of trees otherwise.
}

\seealso{
\code{\link{as.treeshape}} for more details about the convertion of non-binary trees.
}

\references{ 
\url{http://www.ebi.ac.uk/goldman-srv/pandit/} for more details about the Pandit database.
  
\url{http://www.treebase.org/treebase/} for informations about the TreeBASE database. 
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.f}>\cr
Olivier François <\email{olivier.francois@imag.fr}>
}

\examples{
if(!inherits(try(open(url("http://www.google.com")), silent = TRUE),
   "try-error")) {
## Sackin's index of a tree within TreeBASE.
sackin(dbtrees(db="treebase", tree=715, class="treeshape"))

## Colless' index of a tree within Pandit
colless(dbtrees(db="pandit", tree=1, class="treeshape"))
  
## Collects a tree without printings: 
plot(dbtrees(db="pandit", tree=709, class="treeshape", quiet=TRUE))

## Collects a list of trees :
trees=dbtrees(db="pandit", tree=1:5, class="treeshape")
summary(trees[[2]])
  
## Collects a non-binary tree
phy=dbtrees(db="pandit", tree=1, class="phylo")
plot(phy)
tree=dbtrees(db="treebase", tree=741, class="treeshape")
tree=dbtrees(db="treebase", tree=741, class="treeshape", model="yule")
plot(tree)
}
}
\keyword{manip}
\keyword{IO}
