"simulate.brownian" <-
function(phylogeny){
#Add a vector of changes simulated under a Brownian motion model
rands<-rnorm(length(phylogeny$edge.length))
phylogeny$changes<-rands*sqrt(phylogeny$edge.length)
return(phylogeny)
}

