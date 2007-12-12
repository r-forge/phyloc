\name{my.stud.res}
\alias{my.stud.res}
\alias{my.std.res}
\title{ Compute studentized and standardised residuals from pic.lm models }
\description{
Independent contrasts can have heteroscedasticity for several reasons, with the result that 
some points in regression can have undue influence.  These functions are useful for 
identifying scuh points, and are called by \code{pic.lm} and associated functions as part 
of their optional automatic removal of highly influential points.}
\usage{
my.stud.res(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ An object of class \code{"lm"}, hopefully produced by \code{pic.lm}. }
}
\details{
These functions were written because of apparently strange behaviour of the MASS function \code{stdres} 
and hence \code{studres}. }
\value{
numeric vector of either standardised or studentized residuals}
\references{ Venables and Ripley, MASS, pp. 204--205. }

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\seealso{ \code{\link{pic.lm}}, when check.robust = TRUE }


## The function is currently defined as
function(model){
#Written because of strange behaviour of MASS function stdres (and hence studres).
#See Venables and Ripley MASS book, pp 204-5

	std.res<-my.std.res(model)
	p<-length(coef(model))
	n<-length(model$residuals)
	stud.res<-std.res/sqrt((n-p-std.res^2)/(n-p-1))
	
	stud.res
  }
}
\keyword{ manip}
\keyword{ regression}