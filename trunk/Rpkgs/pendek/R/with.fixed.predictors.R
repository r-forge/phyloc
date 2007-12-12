"with.fixed.predictors" <-
function(phylogeny, dataset, response, fixed, predictors, check.robust=FALSE, cutoff=3, use.robust=FALSE, kappa=NULL, filter=NULL){

models<-as.list(predictors)

for(i in 1:length(predictors))
{
	print(paste("Computing model for",predictors[i]))
	models[[i]]<-summary.pic.lm(pic.lm(phylogeny, dataset, response, c(predictors[i],fixed),check.robust=check.robust, cutoff=cutoff, use.robust=use.robust, kappa=kappa))
}

models
}

