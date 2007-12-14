getLambda.internal <-function(zmat, rootnode, rup=NULL, para=.01, eps, combined =TRUE)
{
  int <- zmat[zmat[,2] > rootnode, ]
  term <- zmat[zmat[,2] < rootnode, ]
  nint<-nrow(int)
  nterm<-nrow(term)
  
  betaF <- function(r, t1)
  {
    xf <- (exp(r*t1)-1)/(exp(r*t1)-eps)
    xf;
  }
  Lfunc_tax <- function(p) #likelihood of tax data
  {
    r<-p
   (sum(log(1 - betaF(r, term[1:nterm, 4]))) +
    sum((term[1:nterm, 5]-1)*log(betaF(r, term[1:nterm, 4]))))
  }
  Lfunc_phy <- function(p) #likelihood of phylogenetic data only
  {
    r<-p
    (nint*log(r) - r*sum(int[1:nint,4])
      - sum(log(1-(eps*exp(-r*int[1:nint,3])))))
  }
  Lfunc_comb <- function(p) #combined likelihood
  {
    r<-p
    (sum(log(1 - betaF(r, term[1:nterm, 4]))) +
    sum((term[1:nterm, 5]-1)*log(betaF(r, term[1:nterm, 4])))  #end taxblock#                                            
      + nint*log(r) - r*sum(int[1:nint,4])            #start phyblock#
      - sum(log(1-(eps*exp(-r*int[1:nint,3])))))    
  }
  if(is.null(rup)){
  	rup <- log(1e6)/max(zmat[,3]) #upper bound
  }
  res<-list()
  
  #if/else below necessary to deal with case where 'node' falls at mrca of 
  # terminal branches, thus giving NULL for int, causing the combined 
  # likelihood to fail.
  if (combined ==TRUE){
  	if (nrow(int) == 0)
    	tempres <- optimize(Lfunc_tax, interval=c(0.00001, rup), maximum=TRUE)
 	else
    	tempres <- optimize(Lfunc_comb, interval=c(0.00001, rup), maximum=TRUE)
  }else{
  	 	tempres <- optimize(Lfunc_tax, interval=c(0.00001, rup), maximum=TRUE);
  }
 
 
  res$LH <- tempres$objective
  res$lambda <- tempres$maximum/(1-eps)
  res$r <- tempres$maximum                     
  res$eps <- eps
  res<-as.data.frame(res)
  return(res)
}

