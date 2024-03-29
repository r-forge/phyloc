fitNDR_1rate <- function(phy, eps=0, combined=TRUE)
{
	z <- splitEdgeMatrix(phy, phy$Nnode);
	r1 <- getLambda.internal(z, rootnode = (length(phy$tip.label) +1), rup=NULL, para=.01, eps, combined = combined);
	res<-list();
    res$LH <- r1$LH;
    res$aic <- (-2*r1$LH) + 2;
    res$r <- r1$r;
    res$lambda <- r1$lambda;
    res$eps <-eps;
    res <-as.data.frame(res);
    return(res);
}

