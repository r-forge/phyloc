\encoding{latin1}
\name{gridrowcol}
\alias{gridrowcol}
\title{Complete regular grid analysis}
\description{
This function defines objects to analyse data sets associated with complete regular grid.
}
\usage{
gridrowcol(nrow, ncol, cell.names = NULL)
}
\arguments{
  \item{nrow}{size of the grid (number of rows)}
  \item{ncol}{size of the grid (number of columns)}
  \item{cell.names}{grid cell labels}
}
\value{
Returns a list containing the following items :
    \item{xy}{: a data frame with grid cell coordinates}
    \item{area}{: a data frame with three variables to display grid cells as areas}
    \item{neig}{: an object of class \code{'neig'} corresponding to a neighbouring graph of the grid (rook case)}
    \item{orthobasis}{: an object of class \code{'orthobasis'} corresponding to the analytical solution for the neighbouring graph}
}
\references{
M�ot, A., Chessel, D. and Sabatier, D. (1993) Op�rateurs de voisinage et analyse des donn�es spatio-temporelles.
\emph{in} J.D. Lebreton and B. Asselain, editors. Biom�trie et environnement. Masson, 45-72.

Cornillon, P.A. (1998) \emph{Prise en compte de proximit�s en analyse factorielle et comparative}. Th�se, Ecole Nationale Sup�rieure Agronomique, Montpellier.
}
\author{S�bastien Ollier \email{ollier@biomserv.univ-lyon1.fr} \cr
Daniel Chessel 
}
\seealso{\code{\link{orthobasis}}, \code{\link{orthogram}}, \code{\link{mld}}}
\examples{
w <- gridrowcol(8,5)
par(mfrow = c(1,2))
area.plot(w$area,center = w$xy, graph = w$neig, clab = 0.75)
area.plot(w$area,center = w$xy, graph = w$neig, clab = 0.75,
 label = as.character(1:40))
par(mfrow = c(1,1))

par(mfrow = c(5,8))
for(k in 1:39)
    s.value(w$xy, w$orthobasis[,k], csi = 3, cleg = 0, csub = 2,
     sub = as.character(signif(attr(w$orthobasis, "values")[k],3)),
      incl = FALSE, addax = FALSE, cgr = 0, ylim = c(0,10))
par(mfrow = c(1,1))
}
\keyword{spatial}
