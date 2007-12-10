maxlik.betasplit<-function(phylo,up=10,remove.outgroup=FALSE,confidence.interval="none",prob.conf.inter=.95,size.bootstrap=100)
{
	
	vrais.aldous.fin<-function(i,n,b)
	{
		return(beta(b+i+1,b+n-i+1)/beta(i+1,n-i+1))
	}

	bbalance<-function(phylo)
	{
		return (t(apply(balance(phylo),FUN=function(x){c(x[1],x[1]+x[2])},MARGIN=1)))
	}

	renorm.aldous<-function(n,beta)
	{
		return(sum(sapply(1:(n-1),FUN=vrais.aldous.fin,n=n,b=beta)))
	}

	vrais.aldous.renorm<-function(i,n,beta)
	{
		return(vrais.aldous.fin(i,n,beta)/renorm.aldous(n,beta))
	}

	logvrais.aldous.phylo<-function(b,phylo,remove.outgroup=TRUE,prbootstrap=FALSE)
	{
		if(prbootstrap)
		{
			bal<-bbalance(as.phylo(phylo))
		}
		else
		{
			if(class(phylo)=="treeshape")
				bal<-bbalance(as.phylo(phylo))
	
			if(class(phylo)=="phylo")
				bal<-bbalance(phylo)
	
			if(remove.outgroup)
			{
			if((bal[1,1]<=2)  || ((bal[1,2]-bal[1,1])<=2)) bal<-bal[-1,]
			}
		}
		return (sum(log(apply(bal,FUN=function(x,b){return (vrais.aldous.renorm(x[1],x[2],b))},b=b,MARGIN=1))))
	}

	if(class(phylo)=="treeshape")
		bal<-bbalance(as.phylo(phylo))
	if(class(phylo)=="phylo")
		bal<-bbalance(phylo)
	if(class(phylo)!="phylo" && class(phylo)!="treeshape")
	{
		print("The phylogeny shall be of class phylo or treeshape")
		return
	}
	
	if(remove.outgroup)
	{
		if((bal[1,1]<=2)  || ((bal[1,2]-bal[1,1])<=2)) bal<-bal[-1,]
	}

	nb.tip<-max(bal[1,])
	
	optim_lik_aldous<-function(phylo,remove.outgroup,prbootstrap)
	{
		optimize(f=function(x){logvrais.aldous.phylo(x,phylo,remove.outgroup,prbootstrap)},lower=-2,upper=up,maximum=TRUE)
		#optim(par=0,fn=function(x){logvrais.aldous.phylo(x,phylo,remove.outgroup,prbootstrap)}, 
#method="L-BFGS-B",lower=-1.99,upper=up, hessian = TRUE,control=list(fnscale=-1))
	}
	res<-optim_lik_aldous(phylo,remove.outgroup,prbootstrap=FALSE)

	if (confidence.interval=="bootstrap_param")
	{
		function_aux<-function(n,i)
		{
			if (i == 0 | i == n)
            			return (0)
			else
				return(vrais.aldous.fin(i,n,res$maximum))
				#return(vrais.aldous.fin(i,n,res$par))
		}
		print(res$maximum)
		#print(res$par)
		tree.boot<-rtreeshape(n=size.bootstrap,tip.number=nb.tip,FUN=function_aux)
		thebeta<-sapply(tree.boot,FUN=function(x){optim_lik_aldous(x,remove.outgroup,prbootstrap=TRUE)$maximum})
		#thebeta<-sapply(tree.boot,FUN=function(x){optim_lik_aldous(x,remove.outgroup,prbootstrap=T)$par})
		up.conf<-1-((1-prob.conf.inter)/2)
		low.conf<-(1-prob.conf.inter)/2
		conf_interval<-quantile(thebeta,c(low.conf,up.conf))
	}
	if (confidence.interval=="chi_square")
	{
		x <- list(label=c("beta"),est=c(res$maximum),low=c(-2),upp=c(up))
		conf_interval<-plkhci(x,nlogf=function(x){-logvrais.aldous.phylo(x,phylo,remove.outgroup=FALSE,prbootstrap=FALSE)},label="beta",prob=prob.conf.inter)
		#conf_interval<-res$par+qnorm(c(.025,.975))*sqrt(-res$hessian)
	}
	if (confidence.interval=="none")
	{
		conf_interval<-NULL
	}
	return (list(max_lik=res$maximum,conf_interval=conf_interval))
	#return (list(max_lik=res$par,conf_interval=conf_interval))
}

