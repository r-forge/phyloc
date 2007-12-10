\name{read.inp.data}
\alias{read.inp.data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Read Inp Data Files}
\description{
  Reads one or more inp data files, such as used by PDTREE, of
  the PDAP program bundle, and
  returns an R data frame. Allows to combine several inp files
  and to change the name of variables.
}
\usage{
read.inp.data(input.inp.files, variable.names=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{input.inp.files}{the name (with path if necessary), of
   the inp file(s).}
 \item{variable.names}{an optional vector with the new names for the variables. }
}
\value{
  A data frame (with class pdi.file and data frame) with first column
  the names of tips and remaining columns the data columns from
  the inp file(s).
}
\references{
Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
}

\author{ Diaz-Uriarte, R., and Garland, T., Jr.}

\seealso{
  \code{\link{read.pdi.data}}, \code{\link{read.phylip.data}}, \code{\link{read.sim.data}},  
}

\examples{

%#ifdef unix
# This works under both Unix and Windows.
# First need to find out where the ''Examples'' directory is located.
path.to.example <- paste(.path.package(package="PHYLOGR"),"Examples/",sep="/") 
%#endif


%#ifdef windows
% path.to.example <-
%paste(.path.package(package="phylog"),"Examples\\",sep="\\") # where the example data sets are
%#endif

# a simple case
p49a <- paste(path.to.example,"49lbr.inp",sep="")
data.49a <- read.inp.data(p49a)
data.49a

# two files and rename columns
p49b <- paste(path.to.example,"49hmt.inp",sep="")
data.49.2 <- read.inp.data(c(p49a,p49b),variable.names=c("y","x1","x2","x3"))
data.49.2

# You could jump directly to the call to the function if you
# are willing to enter the path explicitly.
# For example in some Linux systems the following works
# read.inp.data("/usr/lib/R/library/PHYLOGR/Examples/49lbr.inp")
# In Windows, maybe do:
# read.inp.data("c:\\\\progra~1\\\\rw1001\\\\library\\\\PHYLOGR\\\\Examples\\\\49lbr.inp")


}

\keyword{manip}
\keyword{file}
