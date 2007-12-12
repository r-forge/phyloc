`fitDiscrete` <-
function(phy, data, model="ER", data.names=NULL, lambda=FALSE, delta=FALSE, kappa=FALSE, linearchange=FALSE, exponentialchange=FALSE, tworate=FALSE, start.rate=0.01)
{
	td<-treedata(phy, data, data.names, sort=T)

	res<-list()

	for(i in 1:ncol(td$data))
	{
		
		if(lambda + delta+ linearchange + exponentialchange + tworate >1) {
			cat("Cannot handle more than one of (lambda, delta, linearchange, exponentialchange, tworate) at the same time\n");
			return()
		}
	
		if(lambda + delta+ linearchange + exponentialchange + tworate == 0) {
			
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x), model)
			}
						
			nRateCats<-getRateCats(td$data[,i], model)
			
			outTries<-list()
			totalbl<-sum(td$phy$edge.length)
			min=log(0.01/totalbl)
			max=log(10000/totalbl)
			ltry<-numeric(10)
			lsol<-matrix(nrow=10, ncol=nRateCats)
			sp<-numeric(nRateCats)

			for(j in 1:10) {
				sp<-runif(nRateCats, min, max)
				outTries[[j]]<-optim(f, par=sp, method="L",  lower=rep(min, nRateCats), upper=rep(max, nRateCats))
				ltry[j]<-outTries[[j]]$value
				lsol[j,]<-exp(outTries[[j]]$par)

			}
			
			ltd<-ltry-min(ltry)
			gc<-sum(ltd<0.1)				
			b<-min(which(ltry==min(ltry)))
			out<-outTries[[b[1]]]	
			if(gc>1) {out$message="Warning: likelihood surface is flat."}
			
			if(out$convergence!=0) {out$message="Warning: may not have converged to a proper solution."}

			res[[i]]<-list(lnl=-out$value, q=getQ(exp(out$value), nRateCats, model), message=out$message)
			if(!is.null(colnames(td$data))) names(res)[i]<-colnames(td$data)[i] else names(res)[i]<-paste("Trait", i, sep="")
		}
	

	
		if(lambda) {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[1]), lambda=exp(x[2]))
			}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), lambda=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)

		}
		if(delta) {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[1]), delta=exp(x[2]))
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), delta=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)

		}
		if(kappa) {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[1]), kappa=exp(x[2]))
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), kappa=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)

		}
		
		if(linearchange) {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[1]), endRate=exp(x[2]), linear=T)
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), endRate.linear=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)
			
		}
		
		if(exponentialchange) {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[1]), endRate=exp(x[2]))
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), endRate.exponential=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)
			
		}
		
		if(tworate) {
			f<-function(x) {
				likelihoodDiscrete(td$phy, td$data[,i], exp(x[1]), breakPoint=x[2], endRate=exp(x[3]))
				}
			out<-nlm(f, p=c(log(start.rate), 0.5, log(start.rate)))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), breakPoint=out$estimate[2], endRate.tworate=exp(out$estimate[3]), gradient=out$gradient, code=out$code, iterations=out$iterations)
			
		}
	}
	return(res)

}


###Felsenstein's pruning algorithm
likelihoodDiscrete<-function(phy, tip.data, q, model="ER", delta=1, lambda=1,  kappa=1, endRate=1, linear=F, breakPoint=0, f=1, rtt.rescale=0, total.rescale=F, returnFull=F)
{
	
	if(!is.factor(tip.data)) tip.data<-factor(tip.data)
	Q<-getQ(q, n=nlevels(tip.data), model=model)
	if (class(phy) != "phylo")
		stop("object \"phy\" is not of class \"phylo\"");
	#new2old.phylo(phy)->phy ##converts back to old ape tree format with negative values denoting internal nodes
	tmp <- as.numeric(phy$edge)
	nb.tip <- max(tmp) #number of tips
	nb.node <- -min(tmp) #number of internal nodes
	nb.states <- nlevels(tip.data) #numbers of states in character
	l <- matrix(0, nrow=nrow(phy$edge), ncol=nb.states) #makes matrix that will store likelihoods of each character state at each node
	root <- numeric(nb.states) #makes vector that will store root likelihoods
	m <- match(phy$tip.label, names(tip.data)) ##identifies elements of tip.data matrix that corresponds with the tip.label
	if (delta != 1)
		deltaTree(phy, delta) -> phy;
	if (lambda != 1)
		lambdaTree(phy, lambda) -> phy;
	if (kappa != 1)
		kappaTree(phy, kappa) -> phy;
	if (endRate != 1) {
		if(breakPoint!=0) {
			tworateTree(phy, breakPoint, endRate) -> phy;
		} else if(linear==T) {
			linearchangeTree(phy, endRate) -> phy;
		} else exponentialchangeTree(phy, endRate)->phy;
	}


	
	#When comparing deltas across different qs, it might be useful to rescale the total tree length to one	
	if(rtt.rescale!=0)	
		rescaleTree(phy, rtt.rescale) -> phy
		
	new2old.phylo(phy)->phy
	for(i in 1:nrow(phy$edge)) #for each edge
		if(as.numeric(phy$edge[i,2])>0) l[i,tip.data[m[as.numeric(phy$edge[i,2])]]] <- 1.0 #if the edge is connected to a terminal taxon, you set the likelihood of the tip value equal to 1 and all others equal to zero.
		times <- branching.times(old2new.phylo(phy)) #get node to tip distances
		-1*(1:(max(as.numeric(names(times)))-min(as.numeric(names(times)))+1))->names(times)
		times = max(times) - times #convert into root to node tips
		if(total.rescale) {
			sum(phy$edge.length) -> total.tree
			phy$edge.length <- phy$edge.length/total.tree
		}
	while(1) {
		
		if(sum(as.numeric(phy$edge[,2])>0)==2) break 
	
		#obtain ancestors of current tips
		x <- match(phy$edge[,1], phy$edge[,2])[as.numeric(phy$edge[,2])>0] #finds nodes connected to terminal taxa
		#find last node with two tip descendent
		a <- max(x[duplicated(x)])
		t <- which(phy$edge[,1]==phy$edge[a,2])
		bl <- phy$edge.length[t]
		age = times[which(names(times)==phy$edge[a,2])]
		l[a,] <- frag.like(l[t,], bl, Q)

		#next line effectively prunes out the tips just used
		phy$edge[a,2]<-1
		phy$edge[t,2]<-0

	
	}
	t <- which(as.numeric(phy$edge[,2])>0)
	bl <- phy$edge.length[t]
	root <- frag.like(l[t,], bl, Q)
	neglnl=-log(sum(root/nb.states))
	if(returnFull==F) {
		return(neglnl)
	} else return(list(neglnl=neglnl, root=root, l=l))
}

getQ<-function(q, n, model)
{
	if(model=="ER") Q=evenQ(n)*q
	if(model=="SYM") {
		if(length(q)!=n*(n-1)/2) stop("You must supply the correct number of rate categories.")
		Q<-diag(n)
		xx=1
		for(i in 2:n) {
			for(j in 1:(i-1)) {
				Q[i,j]<-Q[j,i]<-q[xx]
				xx<-xx+1
			}
		}	
		for(i in 1:n) diag(Q)[i]<- -sum(Q[i,-i])
	}
	
	if(model=="ARD") {
		if(length(q)!=n*(n-1)) stop("You must supply the correct number of rate categories.")
		Q<-diag(n)
		xx=1
		for(i in 1:n) {
			for(j in (1:n)[-i]) {
				Q[i,j]<-q[xx]
				xx<-xx+1
			}
		}	
		for(i in 1:n) diag(Q)[i]<- -sum(Q[i,-i])
	}
	
	return(Q)
}

getRateCats<-function(data, model)
{
	if(model=="ER") return(1)
	n<-nlevels(factor(data))
	if(model=="SYM") return(n*(n-1)/2)
	if(model=="ARD") return(n*(n-1))
	
}

##evenQ is an internal function of GEIGER
##This function makes a template for the calculate of a rate matrix that will set all transitions to the same value (based on the total number of states).  One can then multiply the resulting matrix by the overall rate to get the rate matrix for a particular analysis.


`evenQ` <-
function(n)

{

	q<--diag(n)

	q[q==0]<-1/(n-1)

	return(q)

}

#This function calculates the likelihood on one branch of a tree

`frag.like` <-
function(tip.like, bl, q)

{

	nb.states<-ncol(tip.like)

	r<-rep(1, nb.states)

	d<-length(bl)

	p<-list(d)

	for(i in 1:d)

		p[[i]]<-MatrixExp.eig(q*bl[i])

	for(i in 1:nb.states)

		for(j in 1:d) 

			r[i]<-r[i]*sum(p[[j]][i,]*tip.like[j,])

	return(r)

}

#This function is required so that the matrix exponentiation never blows up during the likelihood calculation - but it only works for  symmetric matrices (ie evenQ matrices)


`MatrixExp.simple` <-
function(Q)
{
	n<-nrow(Q)
	res<-matrix(0, nrow=n, ncol=n)
	q<-Q[1,2]
	for(i in 1:n)
		res[i, i]<-1/n+(n-1)/n*exp(-n*q)
	res[res==0]<-1/n-1/n*exp(-n*q)
	return(res)
}

## Code lifted from ape function "ace"
MatrixExp.eig<-
function(Q)
{
	 tmp <- eigen(Q, symmetric = FALSE)
	 P1 <- tmp$vectors %*% diag(exp(tmp$values)) %*% solve(tmp$vectors)
	 return(P1)
	
}

