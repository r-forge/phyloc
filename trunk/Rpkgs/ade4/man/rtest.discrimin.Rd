\name{rtest.discrimin}
\alias{rtest.discrimin}
\title{
Monte-Carlo Test on a Discriminant Analysis (in R).
}
\description{
Test of the sum of a discriminant analysis eigenvalues (divided by the rank). Non parametric version of the Pillai's test. It authorizes any weighting.
}
\usage{
\method{rtest}{discrimin}(xtest, nrepet = 99, \dots)
}
\arguments{
  \item{xtest}{an object of class \code{discrimin}}
  \item{nrepet}{the number of permutations}
  \item{\dots}{further arguments passed to or from other methods}
}
\value{
returns a list of class \code{rtest}
}
\author{Daniel Chessel }
\examples{
data(meaudret)
pca1 <- dudi.pca(meaudret$mil, scan = FALSE, nf = 3)
rand1 <- rtest(discrimin(pca1, meaudret$plan$dat, scan = FALSE), 99)
rand1
#Monte-Carlo test
#Observation: 0.3035 
#Call: as.rtest(sim = sim, obs = obs)
#Based on 999 replicates
#Simulated p-value: 0.001 
plot(rand1, main = "Monte-Carlo test")
summary.manova(manova(as.matrix(meaudret$mil)~meaudret$plan$dat), "Pillai")
#                   Df Pillai approx F num Df den Df  Pr(>F)    
# meaudret$plan$dat  3   2.73    11.30     27     30 1.6e-09 ***
# Residuals         16                                          
# ---
# Signif. codes:  0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1 
# 2.731/9 = 0.3034
}
\keyword{multivariate}
\keyword{nonparametric}
