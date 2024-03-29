"gearymoran" <- function (bilis, X, nrepet=999,alter=c("greater", "less", "two-sided")) {
    # bilis doit �tre une matrice
    bilis <- as.matrix(bilis)
    nobs <- ncol(bilis)
    # bilis doit �tre carr�e
    if (nrow(bilis) != nobs) stop ("'bilis' is not squared")
    # bilis doit �tre sym�trique
    bilis <- (bilis + t(bilis))/2
    # bilis doit �tre � termes positifs (voisinages)
    if (any(bilis<0)) stop ("term <0 found in 'bilis'")
    test.names <- names(X)
    X <- data.matrix(X)
    if (nrow(X) != nobs) stop ("non convenient dimension")
    nvar <- ncol(X)
    res <- .C("gearymoran",
        param = as.integer(c(nobs,nvar,nrepet)),
        data = as.double(X),
        bilis = as.double(bilis),
        obs = double(nvar),
        result = double (nrepet*nvar),
        obstot = double(1),
        restot = double (nrepet),
        PACKAGE="ade4"
    )
    res <- as.krandtest(obs=res$obs,sim=matrix(res$result,ncol=nvar, byr=TRUE),names=test.names,alter=alter)
    return(res)       
}
