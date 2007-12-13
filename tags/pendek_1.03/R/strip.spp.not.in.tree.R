"strip.spp.not.in.tree" <-
function(dataset, phylogeny){

#Written by Andy Purvis
#Omits species from a data frame (dataset) that are not present in a phylogeny
#Presumes (as do other functions here) that species names are in the first column of dataset

#requires ape
old<-length(dataset[,1])
matches<-match(dataset[,1],phylogeny$tip.label,nomatch=0)
dataset<-subset(dataset,matches>0)
new<-length(dataset[,1])
print(paste("Dataset stripped from",old,"to",new,"species."))

dataset
}

