"abst" <-
function(kappa, phylogeny, dataset, variable){

#Written by Andy Purvis
#Finds absolute value of t for regression of absolute scaled contrasts on scaling factor (SD)

species.data<-dataset[,names(dataset)==variable]
names(species.data)<-dataset[,1]
phylogeny$edge.length<-phylogeny$edge.length^kappa
dx<-pic(species.data,phylogeny,var.contrasts=TRUE)
SD<-sqrt(dx[,2])
model1<-lm(abs(dx[,1])~SD)
T<-summary(model1)$coefficient[2,3]
abst<-abs(T)

abst
}

