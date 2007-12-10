\name{Helper.functions}

\alias{num.sim.tips}
\alias{number.of.simulations}
\alias{number.of.tips.inp}
\alias{number.of.tips.pdi}
\alias{number.of.tips.sim}
\alias{scan.inp.file}
\alias{scan.pdi.file}
\alias{scan.simulation.file}
\alias{tips.names.inp}
\alias{tips.names.pdi}

  %- Also NEED an `\alias' for EACH other topic documented here.
\title{ Helper functions}
\description{
  Several helper functions used by the PHYLOGR main functions (i.e.,
  these functions are called from other functions).
  These are all one to three lines functions, which are used in lieu of
  calls to read.table, scan, etc. They are of no immediate use for the
  end user, but might be helpful for further programming.
}
\usage{
}
\arguments{
  }
\value{
  Depends on the helper function; here is a summary:
  \item{num.sim.tips}{number of tips and number of simulations of a simulated data set}
  \item{number.of.simulations}{the number of simulations of a simulated data set}
  \item{number.of.tips.inp}{ number of tips in inp data file}
  \item{number.of.tips.pdi}{ditto for pdi}
  \item{number.of.tips.sim}{ ditto for sim}
  \item{scan.inp.file}{the two columns with data from inp file}
  \item{scan.pdi.file}{ditto for pdi}
  \item{scan.simulation.file}{ditto for sim file}
  \item{tips.names.inp}{the names of tips from an inp file}
  \item{tips.names.pdi}{ditto for pdi}
}

\references{ 

  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
  
}

\author{Ramon Diaz-Uriarte and Theodore Garland, Jr.}


\keyword{misc }
