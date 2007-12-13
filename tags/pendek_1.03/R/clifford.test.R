"clifford.test" <-
function(A, B, ew.wrap=FALSE){
    
    A.nm <- deparse(substitute(A))
    B.nm <- deparse(substitute(B))
    
    # do some double checking
    if(! all(dim(A) == dim(B))) stop("Matrices must have same dimensions")
    if(! all(is.na(A) == is.na(B))) {
        warning("Non-identical missing values: using pairwise complete")
        union.na <- is.na(A) | is.na(B)
        A[union.na] <- NA
        B[union.na] <- NA
    }
    
    # alleged sample size
    n <- sum(!is.na(A))
    
    # centre the matrices around 0
    A <- A - mean(A, na.rm=TRUE)
    B <- B - mean(B, na.rm=TRUE)
    
    # get the variance estimates
    sv.A <- var(as.vector(A),na.rm=TRUE)
    sv.B <- var(as.vector(B), na.rm=TRUE)
    
    # get the pearson correlation 
    rpearson <- cor(as.vector(A), as.vector(B), use="p", method="pearson")
    
    # get the autocorrelation estimates from the matrices
    acf.A <- clifford.acf(A, ew.wrap=ew.wrap)
    acf.B <- clifford.acf(B, ew.wrap=ew.wrap)
    
    
    if(all(acf.A$nok != acf.B$nok)) stop("Unexpectedly non-matching N(k)") # paranoia

    ncc <- sum(acf.A$nok*acf.A$acf*acf.B$acf)
    
    var.xy <- ncc/n^2
    if(var.xy <= 0) var.xy = (sv.A*sv.B)/n # see Richardson et al (1989)
    
    var.r <- var.xy/(sv.A*sv.B)
    ess <- 1+1/var.r
    w <- rpearson*((ess-1)^0.5)
    t <- rpearson*(sqrt((ess-2)/(1-rpearson)))
    p <- 2 * pt(-abs(t), ess-2)

    RET <- list(A=list(name=A.nm, mat=A, acf=acf.A$acf, svar=sv.A), 
                B=list(name=B.nm, mat=B, acf=acf.B$acf, svar=sv.B),
                nok=acf.A$nok, rpearson=rpearson, n=n, ncc=ncc, 
                var.xy=var.xy, var.r=var.r, ess=ess, w=w, t=t,p=p)
                
    class(RET) <- "clifford.test"
    return(RET)    
}

