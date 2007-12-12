"test.all.linear" <-
function(phylogeny, dataset, response, predictors, check.robust=FALSE, cutoff=3, use.robust=FALSE, kappa=NULL, filter=NULL){

for(i in 1:length(predictors))
{
	print(paste("Computing model for",predictors[i]))
	test.linear(phylogeny, dataset, c(response,predictors[i]),check.robust=check.robust, cutoff=cutoff, use.robust=use.robust,kappa=kappa,filter=filter)
}

}

