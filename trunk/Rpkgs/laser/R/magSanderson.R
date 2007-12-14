lambda.stem.ci <- function(tb, r, eps=0)
{
  ci<- list()
  betaF <- (exp(r*tb)-1)/(exp(r*tb)-eps)
  ci$upper <- 1 + log(.025)/log(betaF)
  ci$lower <-  1 + log(.975)/log(betaF)
  ci <-as.data.frame(ci)
  ci
}

#lambda.crown
# uses formula from Magallon and Sanderson 2001, crown group (eqn A6).
# validated with their angiosperm data.
lambda.crown.ms01 <- function(n, tb, eps=0) #n <- extant lineages; tb<-time of basal divergence
{
  res<-list()
  x1<- .5*n*(1-eps^2) + 2*eps
  x2 <- .5*(1-eps)*sqrt(n*((n*eps^2)- 8*eps + 2*n*eps + n))
  x <- (log(x1 + x2)-log(2))/tb
  
  res$r <- x
  res$lambda <- x/(1-eps) 
  res
}

#lambda.stem: again, from Magallon&Sanderson.
lambda.stem.ms01 <- function(n, tb, eps=0)
{
  res<-list()
  x<- (1/tb)*log(n*(1-eps)+eps)
  res$r <- x
  res$lambda <- x/(1-eps) 
  res
}

#lambda ML estimator.  Based on Raup 1985, Magallon & Sanderson 2001 etc.
lambda.stem.ml <- function(n, tb, eps=0)
{
  betaF <- function(r, t1)
  {
    xf <- (exp(r*t1)-1)/(exp(r*t1)-eps)
    xf;
  }
  LF <- function(r)
  {
   -(log(1 - betaF(r, tb)) + (n-1)*log(betaF(r, tb)))
  }
  res <- suppressWarnings(nlm(function(p) LF(p[1]), .05))  
  res$lambda <- res$estimate/(1-eps)
  res$r <- res$estimate
  res
}