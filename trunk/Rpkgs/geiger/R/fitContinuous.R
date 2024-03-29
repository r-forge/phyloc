`fitContinuous` <-
function(phy, data, data.names=NULL, model=c("BM", "OU", "lambda", "kappa", "delta", "EB"), bounds=NULL,  meserr=NULL)
{
	
	# sort is T because sub-functions assume data are in
	# this particular order
	
	model<-match.arg(model)
	
	td<-treedata(phy, data, data.names, sort=T)

	ntax=length(td$phy$tip.label)

	if(is.null(meserr)) {
		me=td$data
		me[]=0
		meserr=me	
	} else if(length(meserr)==1) {
		me=td$data
		me[]=meserr
		meserr=me
	} else if(is.vector(meserr)) {
		if(!is.null(names(meserr))) {
			o<-match(rownames(td$data), names(meserr))
			if(length(o)!=ntax) stop("meserr is missing some taxa from the tree")
			meserr<-as.matrix(meserr[o,])
		} else {
			if(length(meserr)!=ntax) stop("No taxon names in meserr, and the number of taxa does not match the tree")
			me<-td$data
			me[]=meserr
			meserr=me
		}
	} else {
		if(!is.null(rownames(meserr))) {
			o<-match(rownames(td$data), rownames(meserr))
			meserr=meserr[o,]
		} else {
			if(sum(dim(meserr)!=dim(td$data))!=0)
				stop("No taxon names in meserr, and the number of taxa does not match the tree")
			print("No names in meserr; assuming that taxa are in the same order as tree")	
		}
	}

	#--------------------------------
    #---    PREPARE DATA LIST     ---
    #--------------------------------
	ds			<- list()
   		ds$tree 		<- td$phy          # TIP data 
    #--------------------------------
    #--- SET MODEL SPECIFICATIONS ---
    #--------------------------------
    cat("Fitting ", model, "model:\n")
    #-----------------------------
    #---  SET PARAMETER BOUNDS ---
    #-----------------------------
    #---- DEFAULT BOUNDS
    bounds.default			 <- matrix(c(0.00000001, 20, 0.0000001,1, 0.000001, 1, 0.00001, 2.999999, 0.0000001, 50, -3, 0), nrow=6, ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("beta", "lambda", "kappa", "delta", "alpha", "a");
    colnames(bounds.default) <- c("min", "max")

 	#---- USER DEFINED PARAMETER BOUNDS
 	if (is.null(bounds)) {
 		bounds <- bounds.default       # USE DEFAULTS
 	}else{
 		if (class(bounds)!="list"){
 			stop("Please specify user defined parameter bounds as a list()")
 		}else{
 			specified   <- !c(is.null(bounds$beta), is.null(bounds$lambda), 
 							  is.null(bounds$kappa), is.null(bounds$delta),  is.null(bounds$alpha), is.null(bounds$a)
 							  )
 			bounds.user <- matrix(c(bounds$beta, bounds$lambda, bounds$kappa, bounds$delta, bounds$alpha, bounds$a), 
 								  nrow=sum(specified), ncol=2, byrow=TRUE
 								  )
 			rownames(bounds.user) <- c("beta", "lambda", "kappa", "delta", "alpha", "a")[specified]
   	 		colnames(bounds.user) <- c("min", "max")
  
   	 		#----  SET FINAL SEARCH BOUNDS
 			bounds <- bounds.default
 			bounds[specified,] <- bounds.user     # Final Bounds
   		} # END if list
   	}  # END user bound if loop
   	#--------------------------------
    #---   APPEND MODEL SETTINGS  ---
    #--------------------------------
  	ds$bounds <- data.frame(t(bounds))

  	ds$model  <- model
  	#--------------------------------
    #---        FIT MODEL         ---
    #--------------------------------
    result<-list()
    for(i in 1:ncol(td$data)) {
    	ds$data=td$data[,i]
    	ds$meserr=meserr[,i]
  		result[[i]]<-fitContinuousModel(ds, print=print)
  		if(!is.null(colnames(td$data))) names(result)[i]<-colnames(td$data)[i] else names(res)[i]<-paste("Trait", i, sep="")

  	}
  	result
}


`fitContinuousModel` <-
function(ds, print=TRUE)
{
	bounds 	<- ds$bounds
	model 	<- ds$model
	n 		<- length(ds$data)

	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
	
	beta.start<-var(ds$data)/max(branching.times(ds$tree))


	out         <- NULL
	
	y			<- ds$data				# TIP data
	tree		<- ds$tree			# Tree
	meserr		<- ds$meserr
	n			<- length(y)
	

	#----------------------------------
	#-----       DEFAULT FIT      -----
	#----------------------------------
	if (model=="BM") {
		k<-2
	
		vcv<-vcv.phylo(tree)

		start=log(beta.start)
		lower=log(bounds[1,"beta"])
		upper=log(bounds[2,"beta"])
		
		foo<-function(x) {
			vv<-exp(x)*vcv
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			-dmvnorm(y, mu, vv, log=T)
		}
		
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par))

	#----------------------------------
	#-----       LAMBDA ONLY      -----
	#----------------------------------
	} else if (model=="lambda"){
		k<-3
		
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","lambda")])
		upper=log(bounds[2,c("beta","lambda")])
		
		
		foo<-function(x) {


			vcv<-vcv.phylo(tree)

			index			<-	matrix(TRUE, n,n)
			diag(index)		<- FALSE
			vcv[index] 	<- vcv[index]*exp(x[2])
			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))

	#----------------------------------
	#-----        KAPPA ONLY      -----
	#----------------------------------
	} else if (model=="kappa"){
		k<-3
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","kappa")])
		upper=log(bounds[2,c("beta","kappa")])
				
		
		foo<-function(x) {

			t<-kappaTree(tree, kappa=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))


	#----------------------------------
	#-----        DELTA ONLY      -----
	#----------------------------------	
	} else if (model=="delta"){
		
		k<-3
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","delta")])
		upper=log(bounds[2,c("beta","delta")])
		
		foo<-function(x) {

			t<-deltaTree(tree, delta=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$value, beta= exp(o$par[1]), delta=exp(o$par[2]))	
		
	#----------------------------------
	#-----        ALPHA ONLY      -----
	#----------------------------------			
	} else if (model=="OU"){
	## modified 12 dec 07 to call ouMatrix(x) instead of vcv.phylo(ouTree(x))

		k<-3
		
		start=log(c(beta.start, 0.5))
		lower=log(bounds[1,c("beta","alpha")])
		upper=log(bounds[2,c("beta","alpha")])
	
		vcvOrig<-vcv.phylo(tree)
		foo<-function(x) {
			vcv <- ouMatrix(vcvOrig, exp(x[2]))
			
			## t<-ouTree(tree, exp(x[2]))
			##vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		
		outTries<-list()
		
		# First one: try near BM solution
		start=c(log(beta.start), -50)
		outTries[[1]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		# Second one: try one with very strong constraints
		tv<-var(y)
		start=log(c(tv*2000, 1000))
		outTries[[2]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
	
		# Third one: try OUCH answer as a starting point
		ot<- ape2ouch(tree, data=y)
	       w.ouch<- hansen.fit(ot$d, ot$node, ot$ancestor, ot$time, interval=c(0.01,10))
		start=log(c(w.ouch$sigma^2, w.ouch$alpha))
		lower<-start-1
		upper<-start+1
		outTries[[3]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
	
		# Try ten random ones
		for(i in 1:10){
			while(1) {

				lower=c(runif(2, min=-20, max=-1))
				upper=lower+runif(2, min=0, max=10)
				start=c(runif(1, min=lower[1], max=upper[1]), runif(1, min=lower[2], max=upper[2]))
				te<-try(outTries[[i+3]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
				if(class(te)!="try-error") break
				}
				
		}
		
		# Try range of alphas
		atry<- -5:4
		stry<- log(tv*2*exp(atry))
		for(i in 1:10){
			while(1) {

				lower=c(-20, -20)
				upper=c(10, 10)
				start=c(stry[i], atry[i])
				te<-try(outTries[[i+13]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
				if(class(te)!="try-error") break
				}
				
		}
		
		
		
		ntries<-23
		ltry<-numeric(ntries)
		lsol<-matrix(nrow= ntries, ncol=2)
		for(j in 1:ntries) {
				ltry[j]<-outTries[[j]]$value
				lsol[j,]<-exp(outTries[[j]]$par)
			}

		ltd<-ltry-min(ltry)
		b<-min(which(ltry==min(ltry)))

		gc<-which(ltd<0.01)
		us<-lsol[gc,1]
		usc<-sum((us-min(us))>0.01)			
		out<-outTries[[b[1]]]	
		if(usc>1) {out$message="Warning: likelihood surface is flat."}
			
		if(out$convergence!=0) {out$message="Warning: may not have converged to a proper solution."}

		results<-list(lnl=-out$value, beta= exp(out$par[1]), alpha=exp(out$par[2]), convergence=out$convergence, message=out$message, k=k)


	#----------------------------------
	#-----        EB ONLY      -----
	#----------------------------------	
	} else if(model=="EB"){

		k<-3
		start=log(c(beta.start, 0.01))
		lower=c(log(bounds[1,"beta"]),bounds[1,"a"])
		upper=c(log(bounds[2,"beta"]),bounds[2,"a"])
		
		foo<-function(x) {
			t<-exponentialchangeTree(tree, a=x[2])

			vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par[1]), a=o$par[2])	}else{
		stop("Parameters  \"lambda, \"kappa\" and \"delta\" can only be fit one at a time currently")
	}
	
	results$aic<-2*k-2*results$lnl
	results$aicc<-2*k*(n-1)/(n-k-2)-2*results$lnl
	results$k<-k
	return(results) 

}



phylogMean<-function(phyvcv, data) 
{
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	
	return(m1 %*% m2)
	
	}
	
ouMatrix <- function(vcvMatrix, alpha) 
{
## follows Hansen 1997; does not assume ultrametricity (AH 12 dec 07)
## vectorized by LJH
  vcvDiag<-diag(vcvMatrix)
  diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
  diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)

  Tij = diagi + diagj - (2 * vcvMatrix)
    
  vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
  return(vcvRescaled) 
}
    	