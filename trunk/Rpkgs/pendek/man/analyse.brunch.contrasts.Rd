\name{analyse.brunch.contrasts}
\alias{analyse.brunch.contrasts}

\title{ Hypothesis-testing on set of contrasts produced by BRUNCH }
\description{
The \code{\link{brunch}} function produces a set of contrasts; this function uses them to test 
whether the two variables specified show correlated evolution, using a range of non-parametric and 
(with warnings) parametric tests.}
\usage{
analyse.brunch.contrasts(contrasts, model.with.x = FALSE)
}

\arguments{
  \item{contrasts}{ list of contrasts, produced by \code{\link{brunch}} }
  \item{model.with.x}{ logical: is the Brownian motion model to be assumed for the x-variable? }
}
\details{
Performs non-parametric tests (sign test, Wilcoxon signed-ranks test) and parametric tests (t-tests on 
both scaled and unscaled contrasts) on the y-contrasts, with a Shapiro-Wilk normality test to help 
assess trustworthiness of the t-tests.  Optionally also performs regressions of both scaled and 
unscaled contrasts, if \code{model.with.x == TRUE}
}
\value{
a named vector of various statistics}
\references{ Purvis, A. and Rambaut, A. 1995 \emph{Comput. Appls. Biosci.} \bold{11}, 247--251, for algorithms. }

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{There is no error-checking}

\seealso{ \code{\link{find.plums}}, \code{\link{brunch.contrast}}, \code{\link{brunch}}}

## The function is currently defined as
function(contrasts,model.with.x=FALSE){

#Does simple analyses on the y contrasts in a set of contrasts produced by brunch
#If model.with.x = TRUE, additionally regresses y contrasts on x contrasts

dy<-unlist(lapply(contrasts,function(x){x$y.contrast}))
vars<-unlist(lapply(contrasts,function(x){if (!is.null(x$y.contrast)) x$variance}))
non.zero.dy<-dy[dy!=0]

bt<-binom.test(sum(non.zero.dy>0,na.rm=TRUE),sum(!is.na(non.zero.dy)))
cat(paste("\nOf ",bt$parameter," non-zero contrasts in ",names(contrasts[[1]])[3],", ",bt$statistic," are positive\n",sep=""))
cat(paste("Binomial test p-value:",format.pval(bt$p.value),"\n\n"))

wt<-wilcox.test(non.zero.dy)
cat(paste("Wilcoxon signed-rank statistic V:",wt$statistic,"\n"))
cat(paste("Wilcoxon p-value:",format.pval(wt$p.value),"\n\n"))

sdy<-dy/sqrt(vars)
tt<-t.test(sdy)
cat(paste("Mean scaled contrast in ",names(contrasts[[1]])[3],": ",format(tt$estimate),"\n", sep=""))
cat(paste("t-test vs mean = 0; t = ",format(tt$statistic),", d.f. = ",tt$parameter,"\n", sep=""))
cat(paste("t-test p-value = ", format.pval(tt$p.value),"\n", sep=""))

st<-shapiro.test(sdy)
cat(paste("(Shapiro-Wilk normality test: W =",format(st$statistic),", p-value =",format.pval(st$p.value),")\n\n"))

utt<-t.test(dy)
cat(paste("Mean unscaled contrast in ",names(contrasts[[1]])[3],": ",format(utt$estimate),"\n", sep=""))
cat(paste("t-test vs mean = 0; t = ",format(utt$statistic),", d.f. = ",utt$parameter,"\n", sep=""))
cat(paste("t-test p-value = ", format.pval(utt$p.value),"\n", sep=""))

ust<-shapiro.test(dy)
cat(paste("(Shapiro-Wilk normality test: W =",format(ust$statistic),", p-value =",format.pval(ust$p.value),")\n\n"))


stats<-c(bt$parameter,bt$statistic,bt$p.value,wt$statistic,wt$p.value,tt$estimate,tt$statistic,tt$parameter,tt$p.value,
	st$statistic,st$p.value,utt$estimate,utt$statistic,utt$parameter,utt$p.value,ust$statistic,ust$p.value)
names(stats)<-c("n.contrasts","n.positive","binomial.p.value","wilcoxon.v","wilcoxon.p.value","mean.scaled.contrast","scaled.t",
	"scaled.t.df","scaled.t.p.value","scaled.shapiro.w","scaled.shapiro.p.value","mean.unscaled.contrast","unscaled.t","unscaled.t.df","unscaled.t.p.value",
	"unscaled.shapiro.w","unscaled.shapiro.p.value")

if (model.with.x==TRUE)
{
	dx<-unlist(lapply(contrasts,function(x){x$x.contrast}))
	sdx<-dx/sqrt(vars)

	sm<-lm(sdy~sdx-1)
	um<-lm(dy~dx-1)

	reg.stats<-c(coef(summary(sm)),coef(summary(um)))
	names(reg.stats)<-c("scaled.b","scaled.se","scaled.regression.t","scaled.regression.p","unscaled.b","unscaled.se",
		"unscaled.regression.t","unscaled.regression.p")
	stats<-c(stats,reg.stats)

	cat(paste("Scaled regression: b = ",format(reg.stats[1]),", t = ", format(reg.stats[3]), ", p = ",format.pval(reg.stats[4]),"\n",sep=""))
	cat(paste("Uncaled regression: b = ",format(reg.stats[5]),", t = ", format(reg.stats[7]), ", p = ",format.pval(reg.stats[8]),"\n\n",sep=""))

	old.par<-par(no.readonly=TRUE)
	par(mfrow=c(1,2))
	plot(sdy~sdx, main=paste("Scaled",names(contrasts[[1]])[3]), sub=names(contrasts[[1]])[4])
	abline(sm)
	plot(dy~dx, main=paste("Uncaled",names(contrasts[[1]])[3]), sub=names(contrasts[[1]])[4])
	abline(um)
	par(old.par)
  }


invisible(stats)
  }
}
\keyword{htest}
\keyword{regression}
\keyword{hplot}