"strip.missing.data" <-
function(dataset, variable.list){

#Written by Andy Purvis
#Omits species from a dataset that have missing values for any of the variables in variable.list
#variable.list is a vector of variable names of interest

old<-length(dataset[,1])
for (i in 1:length(variable.list))
{
	dataset<-subset(dataset,!is.na(dataset[names(dataset)==variable.list[i]]))
}
new<-length(dataset[,1])
print(paste("Dataset stripped from",old,"to",new,"species."))

dataset
}

