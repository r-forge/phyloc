"optimise.kappa" <-
function(phylogeny, dataset, variable){

#Written by Andy Purvis
#variable is the name of the variable of interest
#Ref: T Garland Jr, P H Harvey & A R Ives (1992) Procedures for the analysis of comparative data using
#phylogenetically independent contrasts. Syst. Biol. 41:18-32. (NB - They didn't call it kappa.)

#requires ape

opt.kappa<-optimize(abst,c(0,2),tol=0.01,phylogeny=phylogeny,dataset=dataset,variable=variable)

print(paste("Optimum kappa for",variable,"=",opt.kappa$minimum,"; abs(t) =",opt.kappa$objective))
opt.kappa$minimum
}

