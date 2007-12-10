\name{lm.phylog}
\alias{lm.phylog}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Linear Models from Simulated Data Sets}
\description{
  Fit a linear model to the data from a read.sim.data.
}
\usage{
lm.phylog(formula, data, max.num=0, weights=NULL, exclude.tips=NULL,
          lapply.size=100)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{formula}{ a formula with the same syntax as for any other linear
   model in R (see help for lm: ?lm). }
 \item{data}{a data frame ---with a particular structure--- with the observations; it will often be the
   name of a data frame created using read.sim.data. }
 \item{max.num}{ if different from 0, maximum number of simulations to analyze}
 \item{weights}{ an optional vector of weights to perform weighted least squares; can be a
   column from the data frame or a vector in the parent environment (of
   same length as any column of data before the data is reduced, if
   appropriate, with exclude.tips and max.num arguments).}
 \item{exclude.tips}{ an optional vector giving the names of tips to
   exclude from the analyses.}
 \item{lapply.size}{ a tuning parameter that can affect the speed of
   calculations; see Details}
}
\details{
  This function uses a loop over lapply calls (I got the idea from
  Venables and Ripley (2000), ''S  programming''). By changing the
  lapply.size you can change the size of the block over which lapply is
  used. Changes can make a difference in speed; for instance, in my
  machine from about 1 sec per simulation for a set with 49 species to less than 0.5
  sec. The default value worked well in my machine, but your mileage
  will vary.
}
\value{A list of class phylog.lm with components
  \item{call}{the function call.}
  \item{Fits}{a data frame with the fitted coefficients for that model; the coefficients
    for the non-simulated (''real'') data correspond to sim.counter=0.}
  \item{MarginalTests}{a data frame with the F-tests.  The first column is the overall F-test for
    the model, the rest are the marginal F-tests, respecting the
    marginality principle or hierarchy of terms included.  The
    coefficients for the non-simulated (''real'') data correspond to sim.counter=0.}
}
\references{
  Diaz-Uriarte, R., and Garland, T., Jr., in prep. PHYLOGR:
  an R package for the analysis of comparative data via Monte Carlo
  simulations and generalized least squares approaches.
  
  Searle, S. R. (1987) \emph{Linear models for unbalanced data.} Wiley

  Venables, W. N. and Ripley, B. D. (1999) \emph{Modern applied
    statistics with S-Plus, 3rd ed.} Springer-Verlag.
  
  Venables, W. N. and Ripley, B. D. (2000) \emph{S programming.} Springer-Verlag.

  
}
\author{Ramon Diaz-Uriarte and Theodore Garland, Jr.}
\note{ The marginal F-tests returned are obtained from drop1, and thus
  respect the marginality principle.  For instance,
    if your model is y \eqn{\sim}{~}  x1 + x2*x3 you will see an F for x1 and an F
    for x2:x3 but no F's for x2 or x3.  Discussion can be found,
    for example, in Venables & Ripley, (1999), ch. 6; see also Searle,
    (1987), ch. 6, for the ANCOVA case.

}



\seealso{ \code{\link{summary.phylog.lm}}, \code{\link{plot.phylog.lm}},\code{\link{read.sim.data}}, \code{\link{drop1}}}

\examples{
data(SimulExample)
ex1.lm <- lm.phylog(y ~ x1+diet, weights=x2, max.num=20,
                    exclude.tips=c("La","Ls"), data=SimulExample)
ex1.lm
summary(ex1.lm)
par(mfrow=c(2,2))
plot(ex1.lm)


}
\keyword{models}%-- one or more ...
