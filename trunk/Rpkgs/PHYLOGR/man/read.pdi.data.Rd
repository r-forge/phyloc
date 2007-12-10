\name{read.pdi.data}
\alias{read.pdi.data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Read Pdi Data Files}
\description{
  Reads one or more pdi data files, such as used by PDTREE, of
  the PDAP program bundle, and
  returns an R data frame. Allows to combine several pdi files
  and to change the name of variables.
}
\usage{
read.pdi.data(input.pdi.files, variable.names=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{input.pdi.files}{the name(s), with path if necessary, of
   the pdi file(s).} 
 \item{variable.names}{an optional vector with the new names for the variables. }
}
\value{
  A data frame (with class pdi.file and data frame) with first column
  the names of tips and remaining columns the data columns from
  the pdi file(s).
}
\references{ Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.}

\author{ Diaz-Uriarte, R., and Garland, T., Jr.}

\seealso{
  \code{\link{read.inp.data}}, \code{\link{read.phylip.data}}, \code{\link{read.sim.data}},  
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
p49a <- paste(path.to.example,"49lbr.pdi",sep="")
data.49a <- read.pdi.data(p49a)
data.49a

# two files and rename columns
p49b <- paste(path.to.example,"49hmt.pdi",sep="")
data.49.2 <- read.pdi.data(c(p49a,p49b),variable.names=c("y","x1","x2","x3"))
data.49.2

# You could jump directly to the call to the function if you
# are willing to enter the path explicitly.
# For example in some Linux systems the following works
# read.pdi.data("/usr/lib/R/library/PHYLOGR/Examples/49lbr.pdi")
# In Windows, maybe do:
# read.pdi.data("c:\\\\progra~1\\\\rw1001\\\\library\\\\PHYLOGR\\\\Examples\\\\49lbr.pdi")

}

\keyword{manip}
\keyword{file}
