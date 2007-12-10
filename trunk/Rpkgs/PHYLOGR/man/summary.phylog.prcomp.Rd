\name{summary.phylog.prcomp}
\alias{summary.phylog.prcomp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Summarize a phylog.prcomp object}
\description{
  A 'method' for objects of class phylog.prcomp. Shows the original
  data, and provides p-values and quantiles of the eigenvalues
  based on the simulated data. There is a print 'method'
  for this summary.
}
\usage{
\method{summary}{phylog.prcomp}(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{object}{ an object of class phylog.prcomp such as one returned from a previous
   call to prcomp.phylog. }
 \item{...}{ further arguments passed to or from other methods
   (currently not used).}
}
\details{
  The p-values are calculated in two different ways. The 'component-wise' ones,
  where the p-value for a particular eigenvalue is (number of simulated data sets
  with eigenvalue larger than original (''real'') data + 1) / (number of simulated
  data sets + 1). With this approach, you can find that the p-value for,
  say, the second eigenvalue is smaller than the first, which
  is not sensible. It only makes sense to examine the second
  eigenvalue if the first one is ''significant'', etc. Thus,
  when considering the significance of the second eigenvalue
  we should account for the value of the first. In other words, there is only
  support against the null hypothesis (of no singificant second
  component) if both the first and the second eigenvalue
  from the observed data set are larger than most of the
  simulated data sets. We can account for what happens with
  the first eigenvalue by computing the p-value of the second
  eigenvalue as the number of simulations in which the second
  eigenvalue is larger than the observed, or the
  first simulated eigenvalue is larger than the observed one,
  or both (so that the only cases that count agains the null are those
  where both the first and second simulated eigenvalues are smaller
  than the observed ones). Therefore, with the second set of p-values,
  the p-values are increasing.

  We also provide data for parallel anlysis as in Horn (1965; see also
  Zwick & Velicer 1986 and Lautenschlager 1989) where each eigenvalue
  is compared to the average eigenvalue (for that factor) obtained from
  a simulation.  These can then be used to construct scree plots showing
  both the original and the simulated data.


  

}
\value{ A list (of class summary.phylog.prcomp) with elements
  \item{call}{the call to function prcomp.phylog.}
  \item{original.eigenvalues}{the eigenvalues
    corresponding to the original (''real'') data set.}
  \item{p.value.component}{the component-wise p-value ---see Details}
  \item{p.value.multiple}{the 'multiple-eigenvalue' p-value; see Details}
  \item{quant.eigenvalue}{the quantiles from the simulated
    eigenvalues; linear interpolation is used. Note that
    these quantiles are in the spirit of the  ''component-wise p.values''.}
  \item{num.simul}{the number of simulations used in the analyses}
}


\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.

  Horn, J. L. (1965) A rationale and test for the number of factors in
  factor analysis. \emph{Psychometrica}, 30, 179-185.
  
  Krzanowski, W. J. (1990) \emph{Principles of multivariate analysis}
  Oxford University Press.

  Lautenschlager, G. J. (1989). A comparison of alternatives to
  conducting Monte Carlo analyses for determining parallel analysis
  criteria. \emph{Multivariate Behavioral Research}, 24, 365-395.
  
  Morrison, D. F. (1990) \emph{Multivariate statistcal methods, 3rd ed.} McGraw-Hill.

  Zwick, W. R., and Velicer, W. F. (1986). Comparison of five rules for
  determining the number of components to retain. \emph{Psychological
    Bulletin}, 99, 432-442. 

}
\author{Ramon Diaz-Uriarte and Theodore Garland, Jr.}

\seealso{ \code{\link{read.sim.data}}, \code{\link{prcomp.phylog}}}

\examples{
data(SimulExample)
ex1.prcomp <- prcomp.phylog(SimulExample[,-11]) #the 11th column is a factor
summary(ex1.prcomp)

}
\keyword{multivariate}%-- one or more ...
\keyword{methods}

