"all.one.predictors" <-
function(phylogeny, dataset, response, predictors,check.robust=FALSE, cutoff=3, use.robust=FALSE, kappa=NULL){

models<-as.list(predictors)

for(i in 1:length(predictors))
{
	print(paste("Computing model for",predictors[i]))
	models[[i]]<-summary.pic.lm(pic.lm(phylogeny, dataset, response, predictors[i],check.robust=check.robust, cutoff=cutoff, use.robust=use.robust, kappa=kappa))
}

models
}

