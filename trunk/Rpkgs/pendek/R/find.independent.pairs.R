"find.independent.pairs" <-
function(phylogeny){

#Finds all n phylogenetically independent pairs of taxa in a fully binary phylogeny
#Returns an n x 2 matrix of the tip names

n<-floor(length(phylogeny$tip.label)/2) #Number of comparisons
pairs<-matrix(NA, nrow=n ,ncol=2)
count<-0
while (count<n)
{
	some.pairs<-find.cherries(phylogeny)
	how.many<-length(some.pairs[,1])
	pairs[(count+1):(count+how.many),]<-some.pairs
	count<-count+how.many
	to.drop<-some.pairs
	if (length(phylogeny$tip.label)-1!=length(to.drop)) phylogeny<-drop.tip(phylogeny, to.drop)
}

pairs
}

