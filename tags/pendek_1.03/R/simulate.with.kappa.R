"simulate.with.kappa" <-
function(phylogeny,kappa=1){
#Add a vector of changes simulated under a modified Brownian model in which
#  branch lengths are raised to a power, kappa.
rands<-rnorm(length(phylogeny$edge.length))
phylogeny$changes<-rands*sqrt(phylogeny$edge.length^kappa)
return(phylogeny)
}

