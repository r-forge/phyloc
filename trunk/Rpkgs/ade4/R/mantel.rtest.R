"mantel.rtest" <- function (m1, m2, nrepet = 99) {
    if (!inherits(m1, "dist")) 
        stop("Object of class 'dist' expected")
    if (!inherits(m2, "dist")) 
        stop("Object of class 'dist' expected")
    n <- attr(m1, "Size")
    if (n != attr(m2, "Size")) 
        stop("Non convenient dimension")
    permutedist <- function(m) {
        permutevec <- function(v, perm) return(v[perm])
        m <- as.matrix(m)
        n <- ncol(m)
        w0 <- sample(n)
        mperm <- apply(m, 1, permutevec, perm = w0)
        mperm <- t(mperm)
        mperm <- apply(mperm, 2, permutevec, perm = w0)
        return(as.dist(t(mperm)))
    }
    mantelnoneuclid <- function(m1, m2, nrepet) {
        obs <- cor(unclass(m1), unclass(m2))
        if (nrepet == 0) 
            return(obs)
        perm <- matrix(0, nrow = nrepet, ncol = 1)
        perm <- apply(perm, 1, function(x) cor(unclass(m1), unclass(permutedist(m2))))
        w <- as.rtest(obs = obs, sim = perm, call = match.call())
        return(w)
    }
    if (is.euclid(m1) & is.euclid(m2)) {
        tab1 <- pcoscaled(m1)
        obs <- cor(dist.quant(tab1, 1), m2)
        if (nrepet == 0) 
            return(obs)
        perm <- rep(0, nrepet)
        perm <- unlist(lapply(perm, function(x) cor(dist(tab1[sample(n), 
            ]), m2)))
        w <- as.rtest(obs = obs, sim = perm, call = match.call())
        return(w)
    }
    w <- mantelnoneuclid(m1, m2, nrepet = nrepet)
    return(w)
}
