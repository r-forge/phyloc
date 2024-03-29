"mld"<- function (x, orthobas, level, na.action = c("fail", "mean"), plot=TRUE, dfxy = NULL, phylog = NULL,  ...) 
{

# on fait les v�rifications sur x
if (!is.numeric(x)) 
        stop("x is not numeric")
nobs <- length(x)
if (any(is.na(x))) {
        if (na.action == "fail") 
            stop(" missing values in 'x'")
        else if (na.action == "mean") 
            x[is.na(x)] <- mean(na.omit(x))
        else stop("unknown method for 'na.action'")
    }

# on fait les v�rifications sur orthobas (class, dimension, orthogonalit�, orthonormalit�)
if (!inherits(orthobas, "data.frame")) stop ("'orthobas' is not a data.frame")
    if (nrow(orthobas) != nobs) stop ("non convenient dimensions")
    if (ncol(orthobas) != (nobs-1)) stop (paste("'orthobas' has",ncol(orthobas),"columns, expected:",nobs-1))

vecpro <- as.matrix(orthobas)

w <- t(vecpro/nobs)%*%vecpro
    if (any(abs(diag(w)-1)>1e-07)) {
        stop("'orthobas' is not orthonormal for uniform weighting")
    }
    
diag(w) <- 0
    if ( any( abs(as.numeric(w))>1e-07) ) stop("'orthobas' is not orthogonal for uniform weighting")

# on calcule les diff�rents vecteurs associ�s � la d�composition orthonormale de la variable

    # si x n'est pas centr�e, on la centre pour la pond�ration uniforme
    if (mean(x)!=0)
            x <- x-mean(x)
    
    # on calcul les coefficients de corr�lation entre la variable et les vecteurs de la base
    coeff <- t(vecpro/nobs)%*%as.matrix(x)
    
    # on calcul les vecteurs associ�s � la d�composition et au facteur level
    if (!is.factor(level))
            stop("'level' is not a factor")
            if (length(level) != (nobs-1)) 
                    stop (paste("'level' has",length(level),"values, expected:",nobs-1))
    res <- matrix(0, nrow = nobs, ncol = nlevels(level))
    coeff <- split(coeff, level)
    vecpro <- as.data.frame(t(vecpro))
    vecpro <- split(vecpro, level)
    for (i in 1:nlevels(level)) 
        res[,i] <- t(vecpro[[i]])%*%as.matrix(coeff[[i]])
    res <- as.data.frame(res)
    names(res) <- paste("level", levels(level), sep=" ")

 
# on fait les sorties graphiques si elles sont demand�es: c'est pas parfait mais c'est pour donner une id�e
if (plot==TRUE){
    # rajouter les donn�es circulaires
    if (is.ts(x)){
        # pour les s�ries temporelles
        u <- attributes(x)$tsp
        tab <- ts(res, start = u[1], end = u[2], frequency = u[3])
        tab <- ts.union(x, tab)
        u <- range(tab)
        opar <- par(mfrow = par("mfrow"), mar = par("mar"))
        on.exit(par(opar))
        mfrow <- n2mfrow(nlevels(level)+1)
        par(mfrow = mfrow)
        par(mar = c(2.5, 5, 1.5, 0.6))
        plot.ts(x, ylim = u, ylab = "x", main = "multi-levels decomposition")
        for (i in 1:nlevels(level))
                plot(tab[,i+1], ylim = u, ylab = names(res)[i], main = "")
        }
        
    if (is.vector(x)){
        if (!is.null(dfxy)){
            # pour les donn�es 2 D
            opar <- par(mfrow = par("mfrow"), mar = par("mar"))
            on.exit(par(opar))
            mfrow <- n2mfrow(nlevels(level)+1)
            par(mfrow = mfrow)
            par(mar = c(0.6, 2.6, 0.6, 0.6))
            s.value(dfxy, x, sub = "x", ...)
            for (i in 1:nlevels(level))
                if (max((1:(nobs-1))[level == levels(level)[i]])<(nobs/2)){
                    s.image(dfxy, res[,i])
                    s.value(dfxy, res[,i], sub = names(res)[i], add.plot=TRUE, ...)
                    }
                else
                    s.value(dfxy, res[,i], sub = names(res)[i], ...)
            }
        else {
            if (!is.null(phylog)){
                # pour les donn�es associ�es � une phylog�nie
                tab <- cbind.data.frame(x, res)
                row.names(tab) <- names(phylog$leaves)
                table.phylog(tab, phylog, ...)                
                }
            else {
                # pour les transects
                par(mfrow = c(nlevels(level)+1,1))
                par(mar = c(2, 5, 1.5, 0.6))
                u <- range(cbind(x, res))
                w <- trunc(u)
                w <- c(w[1],0,w[2])
                plot(x, type="h", ylim = u, axes = FALSE, ylab = "x", main = "multi-levels decomposition")
                axis(side = 2, at = w, labels = as.character(w))
                for (i in 1:nlevels(level)){
                    plot(res[,i], type="h", ylim = u, axes = FALSE, ylab = names(res)[i], main = "")
                    axis(side = 2, at = w, labels = as.character(w))
                    }
                v <- seq(0, nobs, by = (nobs/10))
                axis(side=1, at = v, labels = as.character(v))
                }        
            }        
        }
    }
return(res)  
}

#############################################################################
haar2level <- function(x){
# cette fonction calcul le facteur level pour lequel l'analyse mld correspond
# � l'analyse mra de la library(waveslim)

# on v�rifie que x=2**a
a <- log(length(x))/log(2)
b <- floor(a)
if ((a-b)^2>1e-10) stop ("Haar is not a power of 2")

#on construit les J niveaux de d�composition
u <- LETTERS[1:a]
v <- rep(2,a)**(0:(a-1))
level <- rep(u, v)
level <- as.factor(level)
return(level)
}
