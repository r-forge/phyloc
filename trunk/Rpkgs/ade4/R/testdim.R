"testdim" <- function (dudi, ...) UseMethod("testdim")

"testdim.pca" <-
  function(dudi,nrepet=99,nbax=dudi$rank,alpha=0.05,...){
    if (!inherits(dudi, "dudi"))
      stop("Object of class 'dudi' expected")
    if (!inherits(dudi, "pca"))
      stop("Object of class 'pca' expected")
    appel <- as.list(dudi$call)
    appel$scale <- eval(appel$scale,sys.frame(0))
    appel$center <- eval(appel$center,sys.frame(0))
    if (is.null(appel$scale)) appel$scale <- TRUE
    if (is.null(appel$center)) appel$center <- TRUE
    if (!(is.logical(appel$center))) stop("Not implemented for decentred PCA")
    if (!(appel$center == TRUE  && appel$scale == TRUE))
      stop("Only implemented for PCA on correlation matrix (center=TRUE and scale=TRUE)")
    X <- as.matrix(dudi$tab)
    if (!(identical(all.equal(dudi$lw,rep(1/nrow(X), nrow(X))),TRUE)))
      stop("Not implemented for non-uniform row weights")
    if (!(identical(all.equal(dudi$cw,rep(1, ncol(X))),TRUE)))
      stop("Not implemented for non-uniform column weights")
    if (nbax<1)
      stop("Incorrect number of axes")
    nbax <- ifelse(nbax>min(nrow(X),ncol(X)),min(nrow(X),ncol(X)),nbax)
    res <- list()
    res <- .C("testdimRVpca", as.double(t(X)), as.integer(nrow(X)), as.integer(ncol(X)), as.integer(nrepet),nbax=as.integer(nbax),sim=as.double(rep(0,nbax*nrepet)),obs=as.double(rep(0,nbax)),PACKAGE="ade4")[c("obs","sim")]
    res$sim <- matrix(res$sim[1:(nbax*nrepet)],nrepet,nbax,byrow=TRUE)
    res$obs <- res$obs[1:nbax]
    res <- as.krandtest(sim=res$sim,obs=res$obs,names=paste("Axis", 1:length(res$obs)),call=match.call())
    
    nb <- which(res$pvalue>alpha)
    if(length(nb)==0) {res$nb <- length(res$obs)} else {res$nb <- min(nb)-1}
    nb2 <- which(res$pvalue>(alpha/1:length(res$obs)))
    if(length(nb2)==0) {res$nb.cor <- length(res$obs)} else {res$nb.cor <- min(nb2)-1}
    return(res)
  }
