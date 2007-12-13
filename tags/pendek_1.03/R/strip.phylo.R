"strip.phylo" <-
function(phylogeny, dataset){

#Written by Andy Purvis
#Drops tips from a phylogeny that are not present in the data set

#requires ape
old<-sum(phylogeny$edge>0)
not.in.data<-setdiff(phylogeny$tip.label,dataset[,1])
if (length(not.in.data)>0) phylogeny<-drop.tip(phylogeny,not.in.data)
new<-sum(phylogeny$edge>0)
print(paste("Phylogeny stripped from",old,"to",new,"species."))

phylogeny
}

