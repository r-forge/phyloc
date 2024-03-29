"kplot.mfa" <- function (object, xax = 1, yax = 2, mfrow = NULL, which.tab = 1:length(object$blo),
    row.names = FALSE, col.names = TRUE, traject = FALSE, permute.row.col = FALSE, 
    clab = 1, csub = 2, possub = "bottomright", ...) 
{
    if (!inherits(object, "mfa")) 
        stop("Object of type 'mfa' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.tab))
    par(mfrow = mfrow)
    if (length(which.tab) > prod(mfrow)) 
        par(ask = TRUE)
    for (ianal in which.tab) {
        coolig <- object$lisup[object$TL[, 1] == ianal, c(xax, yax)]
        coocol <- object$co[object$TC[, 1] == ianal, c(xax, yax)]
        if (permute.row.col) {
            auxi <- coolig
            coolig <- coocol
            coocol <- auxi
        }
        cl <- clab * row.names
        if (cl > 0) 
            cpoi <- 0
        else cpoi <- 2
        s.label(coolig, clab = cl, cpoi = cpoi)
        if (traject) 
            s.traject(coolig, clab = 0, add.p = TRUE)
        born <- par("usr")
        k1 <- min(coocol[, 1])/born[1]
        k2 <- max(coocol[, 1])/born[2]
        k3 <- min(coocol[, 2])/born[3]
        k4 <- max(coocol[, 2])/born[4]
        k <- c(k1, k2, k3, k4)
        coocol <- 0.7 * coocol/max(k)
        s.arrow(coocol, clab = clab * col.names, add.p = TRUE, 
            sub = object$tab.names[ianal], possub = possub, csub = csub)
    }
}
