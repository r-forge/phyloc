"optimal.pic" <-
function(phylogeny, dataset, variable){

#Written by Andy Purvis
#Returns set of contrasts computed with the optimally transformed branch lengths
#variable is the name of the variable of interest

#requires ape

kappa<-optimise.kappa(phylogeny, dataset, variable)
phylogeny$edge.length<-phylogeny$edge.length^kappa
species.data<-dataset[,names(dataset)==variable]
names(species.data)<-dataset[,1]
dx<-pic(species.data,phylogeny)
dx
}

