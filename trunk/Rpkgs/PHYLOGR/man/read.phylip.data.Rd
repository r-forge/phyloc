\name{read.phylip.data}
\alias{read.phylip.data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Read Phylip Infile Data Files}
\description{
  Reads one file with tip data, such as a PHYLIP infile, and
  returns an R data frame. The files have to follow the PHYLIP
  standard of having a first line with the number of species and traits.
  There can be multiple traits, as allowed in PHYLIP. Data for a species
  can extend over several lines, as for
  PHYLIP's sequential data format for continuous traits. It is assumed
  that all traits are cuantitative.
}
\usage{
read.phylip.data(input.phylip.file, variable.names=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{input.phylip.file}{the name, with path if necessary, of
   the infile.} 
 \item{variable.names}{an optional vector with the new names for the variables. }
}
\value{
  A data frame (with class phylip.file and data frame) with first column
  the names of tips and remaining columns the data columns from
  the phylip file.
}
\references{ Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.}

\note{
  The format of PHYLIP's infiles is not exactly the same as the TIP
  files produced from PDTREE. First, PHYLIP's infiles can contain an
  arbitrary number of traits, whereas PDTREE's TIP files only have
  two. Second, PHYLIP's infiles have a first line with the number of tips
  and the number of traits, separated by blanks. Third, the first field
  or column of PHYLIP's infiles must be ten characters wide; if the tips
  names are shorter, they must be padded with blanks (see PHYLIP's
  documentation). This limitation does not apply to read.phylip.data,
  but you might want to follow it if you plan to use PHYLIP.
  }
  
\author{ Diaz-Uriarte, R., and Garland, T., Jr.}

\seealso{
  \code{\link{read.inp.data}}, \code{\link{read.pdi.data}},  \code{\link{read.sim.data}},  
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

lacertid.data.name <- paste(path.to.example,"LacertidData.PhylipInfile",sep="")
lacertid.data <- read.phylip.data(lacertid.data.name,variable.names=c("svl", "svl.matur",
                                  "hatsvl", "hatweight", "clutch.size",
                                  "age.mat", "cl.freq", "xx"))
lacertid.data

# You could jump directly to the call to the function if you
# are willing to enter the path explicitly.
# For example in some Linux systems the following works
# read.phylip.data("/usr/lib/R/library/PHYLOGR/Examples/LacertidData.PhylipInfile")
# In Windows, maybe do:
# read.pdi.data("c:\\\\progra~1\\\\rw1001\\\\library\\\\PHYLOGR\\\\Examples\\\\LacertidData.PhylipInfile")

}

\keyword{manip}
\keyword{file}
