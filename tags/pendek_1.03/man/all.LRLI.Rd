\name{all.LRLI}
\alias{all.LRLI}
\title{ Produce several variants of LRLI plots on phylogenetically structured data }
\description{
Plots and models log(Rate) vs log(Interval) (LRLI) in three ways: 1. Highly non-independent analysis, 
using all-by-all comparisons (twice!); 2. Independent analysis of cherries only; 3. Independent analysis 
of maximum number of independent pairs.}
\usage{
all.LRLI(phylogeny, data, plot = "n")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ object of class \code{phylo}, assumed to be binary }
  \item{data}{ numeric vector of data for all the terminals in \code{phylogeny}. NOTE: assumed to be
  in same order as \code{phylogeny$tip.label} but this is not checked! A later version will be
  restructured to prevent the obvious scope for terrible mistakes.}
  \item{plot}{ is plot required? Default = n }
}

}
\value{
list with three components, each a regression model of class \code{lm}:
  \item{model1 }{Regression using all tip-to-tip distances - highly non-independent, especially as each 
  distance is actually used twice}
  \item{model2 }{Regression using cherries only}
  \item{model3 }{Regression using the maximum number of independent pairs}
Plots are optionally produced.
}
\references{Gingerich, P. D. 1993 \emph{American Journal of Science} \bold{293A}, 453--478, 
for LRLI plots, though the phylogenetic context is novel. }

\author{Andy Purvis}

\section{Version }{1.3}
\section{Date }{19 May 2006}

\section{Warning }{At present there is no error checking}
\section{Warning }{If species values are not in the same order as \code{phylogeny$tip.label}, the function will
return plausible garbage}


\seealso{ \code{\link{matched.pairs.LRLI}}}


## The function is currently defined as
function(phylogeny, data,plot="n"){

#Produces Gingerich-type log-rate log-interval plot for a phylogeny
#Returns three lm objects as a list
#  1. Highly nonindependent analysis using all-by-all comparisons (twice!)
#  2. Independent analysis of cherries only
#  3. Independent analysis of maximum number of pairs

#10/5/06: Bug fix to prevent crash when two species have exactly the same trait value and
#so have zero rate.  See also matched.pairs.LRLI

#10/5/06: Change to reflect change in matched.pairs.LRLI, which now returns data rather than a model.

#11/5/06: Realisation that the function assumes the species data are in the same order as tip.label.
#I might want to consider either using row.names or passing a data frame instead of a vector.

phylogeny.mat<-cophenetic(phylogeny)
data.mat<-as.matrix(dist(data))
rate.mat<-data.mat/phylogeny.mat

#Non-phylogenetic analysis
#Note that, in addition to the obvious problem of nonindependence, both halves of the matrix are used
#so there are twice as many points as there should be
rate<-as.vector(rate.mat)
rate[rate==0]<-min(rate[rate>0])/10 #Hack
interval<-as.vector(phylogeny.mat)

log10.rate<-log10(rate)
log10.interval<-log10(interval)

model1<-lm(log10.rate~log10.interval)

if (plot == "y")
{
	par(mfrow=c(3,1))
	plot(log10.rate~log10.interval,ylab="Log10(Rate)",xlab="Log10(Interval)",cex.lab=1.3)
	abline(model1)
  }

#Phylogenetic analysis of cherries only
cherries<-find.cherries(phylogeny)
ds<-matched.pairs.LRLI(phylogeny, phylogeny.mat, cherries, data.mat)
model2<-lm(ds$Log10.Rate~ds$Log10.Interval)

if (plot == "y")
{
	plot(ds$Log10.Rate~ds$Log10.Interval,ylab="Log10(Rate)",xlab="Log10(Interval)",cex.lab=1.3)
	abline(model2)
  }

#Phylogenetic analysis of all independent pairs
all.matched.pairs<-find.independent.pairs(phylogeny)
ds<-matched.pairs.LRLI(phylogeny, phylogeny.mat, all.matched.pairs, data.mat)
model3<-lm(ds$Log10.Rate~ds$Log10.Interval)

if (plot == "y")
{
	plot(ds$Log10.Rate~ds$Log10.Interval,ylab="Log10(Rate)",xlab="Log10(Interval)",cex.lab=1.3)
	abline(model3)
	par(mfrow=c(1,1))
  }

to.return<-list(model1,model2,model3)

return(to.return)

  }
}
\keyword{ regression }
\keyword{ plot }
