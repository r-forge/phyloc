`fit.discrete` <-
function(phy, data, data.names=NULL, lambda=FALSE, delta=FALSE, kappa=FALSE, linearchange=FALSE, exponentialchange=FALSE, tworate=FALSE, start.rate=0.01)
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
				likelihood.discrete(td$phy, td$data[,i], exp(x))
			}
			outTries<-list()
			outTries[[1]]<-optim(f, par=log(0.01), method="L",  lower=-50, upper=50)
			outTries[[2]]<-optim(f, par=log(0.1), method="L",  lower=-50, upper=50)
			outTries[[3]]<-optim(f, par=log(1), method="L",  lower=-50, upper=50)
			outTries[[4]]<-optim(f, par=log(10), method="L",  lower=-50, upper=50)

			l4<-numeric(4)
			for(j in 1:4)
				l4[j]<-outTries[[j]]$value					
			b<-min(which(l4==min(l4)))
			out<-outTries[[b]]	
			
			res[[i]]<-list(lnl=-out$value, q=exp(out$par), convergence=out$convergence, message=out$message)
		}
	

	
		if(lambda) {
			f<-function(x) {
				likelihood.discrete(td$phy, td$data[,i], exp(x[1]), lambda=exp(x[2]))
			}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), lambda=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)

		}
		if(delta) {
			f<-function(x) {
				likelihood.discrete(td$phy, td$data[,i], exp(x[1]), delta=exp(x[2]))
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), delta=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)

		}
		if(kappa) {
			f<-function(x) {
				likelihood.discrete(td$phy, td$data[,i], exp(x[1]), kappa=exp(x[2]))
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), kappa=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)

		}
		
		if(linearchange) {
			f<-function(x) {
				likelihood.discrete(td$phy, td$data[,i], exp(x[1]), endRate=exp(x[2]), linear=T)
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), endRate.linear=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)
			
		}
		
		if(exponentialchange) {
			f<-function(x) {
				likelihood.discrete(td$phy, td$data[,i], exp(x[1]), endRate=exp(x[2]))
				}
			out<-nlm(f, p=rep(log(start.rate), 2))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), endRate.exponential=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations)
			
		}
		
		if(tworate) {
			f<-function(x) {
				likelihood.discrete(td$phy, td$data[,i], exp(x[1]), breakPoint=x[2], endRate=exp(x[3]))
				}
			out<-nlm(f, p=c(log(start.rate), 0.5, log(start.rate)))
			res[[i]]<-list(lnl=-out$minimum, q=exp(out$estimate[1]), breakPoint=out$estimate[2], endRate.tworate=exp(out$estimate[3]), gradient=out$gradient, code=out$code, iterations=out$iterations)
			
		}
	}
	return(res)

}



