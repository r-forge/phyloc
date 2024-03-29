\encoding{latin1}
\name{euro123}
\alias{euro123}
\docType{data}
\title{Triangular Data}
\description{
This data set gives the proportions of employement in the primary, secondary and tertiary 
sectors for 12 European countries in 1978, 1986 and 1997. 
}
\usage{data(euro123)}
\format{
    \code{euro123} is a list of 4 components.
    \describe{
    \item{in78}{is a data frame with 12 rows and 3 variables.}
    \item{in86}{: idem in 1986}
    \item{in97}{: idem in 1997}
    \item{plan}{is a data frame with two factors to both organize the 3 tables.}
    }
}
\source{
Encyclopaedia Universalis, Symposium, Les chiffres du Monde. Encyclopaedia Universalis, Paris. 519.\cr\cr
Université de Barcelone : \url{http://www.ub.es/medame/nutstat1.html}
}
\examples{
data(euro123)
par(mfrow = c(2,2))
triangle.plot(euro123$in78, addaxes = TRUE)
triangle.plot(euro123$in86, addaxes = TRUE)
triangle.plot(euro123$in97, addaxes = TRUE)
triangle.biplot(euro123$in78, euro123$in97)
par(mfrow = c(1,1))}
\keyword{datasets}
