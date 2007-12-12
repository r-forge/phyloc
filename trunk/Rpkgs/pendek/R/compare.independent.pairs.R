"compare.independent.pairs" <-
function(phylogeny, dataset, variable.list,cherries.only=FALSE){

#Generates BRUNCH-like contrasts in vars from a fully binary phylogeny.

dataset<-strip.spp.not.in.tree(dataset,phylogeny)
dataset<-strip.missing.data(dataset,variable.list)
phylogeny<-strip.phylo(phylogeny,dataset)
nc<-floor(length(dataset[,1])/2)
if (cherries.only==FALSE) pairs<-find.independent.pairs(phylogeny) else pairs<-find.cherries(phylogeny)
names<-list(NULL,variable.list)
contrasts<-array(data=NA,dim=c(length(pairs)/2,length(variable.list)),dimnames=names)
variances<-array(data=NA,dim=length(pairs)/2)

for (i in 1:length(variable.list))
{
	species.data<-dataset[,names(dataset)==variable.list[i]]
	names(species.data)<-dataset[,1]
	pairs.indices<-matrix(match(pairs,names(species.data)),ncol=2)
	contrasts[,i]<-species.data[pairs.indices[,1]]-species.data[pairs.indices[,2]]
}


phylogeny.mat<-dist.phylo(phylogeny)

for (i in 1:length(contrasts)/2)
{
	variances[i]<-phylogeny.mat[pairs.indices[i,1],pairs.indices[i,2]]
}


to.return<-list(contrasts,variances)
names(to.return)<-c("unscaled.contrasts","variances")
to.return
}

