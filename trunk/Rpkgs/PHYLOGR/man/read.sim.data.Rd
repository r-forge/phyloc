\name{read.sim.data}
\alias{read.sim.data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Read Sim Files From PDSIMUL}
\description{
  Reads ouput file(s) from PDSIMUL ---simulation of phenotypic evolution
  along a phylogeny--- and converts into an R data frame.
  Can add original data from  inp or pdi file(s).
}
\usage{
read.sim.data(sim.files, inp.files = NULL, pdi.files = NULL,
              phylip.file = NULL, variable.names = NULL,
              other.variables = NULL, max.num = 0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{sim.files}{ the sim file(s), with path if needed. }
 \item{inp.files}{ the inp file(s), with path if needed. Might have
   already been processed by read.inp.data. If already read, enter name unquoted.}
 \item{pdi.files}{ the pdi file(s), with path if needed. Might have
   already been processed by read.pdi.data. If already read, enter name
   unquoted. }
 \item{phylip.file}{ the phylip.file, with path if needed. Might have
   already been processed by read.phylip.data. If already read, enter
   name unquoted. }
 \item{variable.names}{ a optional vector of variable names}
 \item{other.variables}{ an optional set of other variables that you want to add to the
   data set. Can be a vector, a matrix, or a data frame. It must be the
   same length as the number of tips in you sim file(s). }
 \item{max.num}{ if different from 0, the number of simulations to process.}
}
\details{
  
  You will almost always want to provide inp or pdi or phylip file(s) since this
  is what will be used for the analyses of simulated data sets. The sim
  and pdi (or inp or phylip) files should match one to one; for example, you might
  have used f1.pdi to obtain f1.sim and f2.pdi to produce f2.sim. Then,
  the order of files ought to be
  
  read.sim.data(c("f1.sim","f2.sim"),pdi.files=c("f1.pdi","f2.pdi"))

  or

  read.sim.data(c("f2.sim","f1.sim"),pdi.files=c("f2.pdi","f1.pdi"))

  but NOT

  read.sim.data(c("f1.sim","f2.sim"),pdi.files=c("f2.pdi","f1.pdi"))

  and NOT

  read.sim.data(c("f2.sim","f1.sim"),pdi.files=c("f1.pdi","f2.pdi"))

  since the last two will yield meaningless results.
  
  Remember  that the number of sim file(s) and pdi (or inp) files must match
  (since that is the only way the number of columns will
  match). 

  This does not apply to PHYLIP infiles, since a PHYLIP infile can
  contain multiple columns.
  
  Inp and Pdi files can not (yet) be mixed in a single call. If you need
  to, you should use read.inp.data and read.pdi.data, and change the
  class of the output data frame.
  
  If you are entering inp files only, you don't need to provide the
  argument name. If you are using pdi files you need to provide the
  pdi.files=.

}
\value{
  A data frame (of class simul.phylog and data.frame) where the first column is called sim.counter, for
  simulation counter (with value 0 for the pdi, or inp, data set),
  second column is called tips, and the rest are data columns (including, if given, the other.variables column).
}
\references{
Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
}
\author{ Ramon Diaz-Uriarte and Thodore Garland, Jr.}

\seealso{ \code{\link{read.pdi.data}}, \code{\link{read.inp.data}}, \code{\link{read.phylip.data}}. There
  are generic functions plot and summary.}


\examples{
%#ifdef unix
# First we need to find where the Examples directory is.
# You could enter it directly (see read.pdi.data for an example).
 path.to.example <- paste(.path.package(package="PHYLOGR"),"Examples/",sep="/") 
%#endif


%#ifdef windows
% path.to.example <-
%paste(.path.package(package="phylog"),"Examples\\",sep="\\") # where the example data sets are
%#endif

# simple example
p49.i <- paste(path.to.example,"49lbr.pdi",sep="") 
p49.s <- paste(path.to.example,"49lbr.sim",sep="")
data.49.s <- read.sim.data(p49.s, pdi.files=p49.i)
data.49.s

# several files, added variables, change column names,
# and limit number of cases

f491s <- paste(path.to.example,"49lbr.sim",sep="")
f492s <- paste(path.to.example,"49hmt.sim",sep="")
f493s <- paste(path.to.example,"49ms.sim",sep="")

f491i <- paste(path.to.example,"49lbr.pdi",sep="")
f492i <- paste(path.to.example,"49hmt.pdi",sep="")
f493i <- paste(path.to.example,"49ms.pdi",sep="")


data.hb <-
      read.sim.data(c(f491s, f492s, f493s), pdi.files=c(f491i, f492i, f493i),
                    variable.names=c("x1","x2","x3","x4","x5","x6"),
                    other.variables=data.frame(
                                     mood=c(rep("good",15),
                                            rep("bad",15),
                                            rep("terrible",19)),
                                     color=c(rep("blue",20),
                                             rep("white",29))),
                     max.num=20)

data.hb
}

\keyword{manip}
\keyword{file}
