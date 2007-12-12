\name{makeNamedVec}
\alias{makeNamedVec}
\title{data.frame column to named vector}
\usage{makeNamedVec(x, colID)}
\arguments{
  \item{x}{a data.frame}
  \item{colID}{column of x to be converted into a vector with names from row.names(x)}
}
\description{
  \code{makeNamedVec} is a utility function for converting a column of a data frame to a vector, keeping the data.frame's row.names as the names of the new vector.
}
\value{
 A vector with names (if x has row.names).
}
\author{Steven Kembel \email{skembel@berkeley.edu}}

\keyword{manip}
