# BAYESIAN TOOLS
# Last update: CH 14.08.2007

# Contents:
# 1. stat.posterior
# 2. MrBayes pipe

# 1: This function plots posterior probabilities as obtained by MrBayes over generations.

stat.posterior <- function(filename, path="/Applications/mrbayes-3.1.2/"){
	filename <- paste(path, filename, sep="")
	B1 <- read.table(filename, skip=2, header=FALSE, fill=TRUE)[,c(1,2)]
	plot(B1$V1, B1$V2, main=paste("Posterior probalilities over generations\n", filename), xlab="Number of generations", ylab="Posterior probabilities", cex=0.8, col="red")
	}


# 2: This function runs MrBayes

mrbayes <- function(alignment, filename, nst=6, rates="invgamma", ngammacat=4, nruns=2, ngen=1000000, printfreq=100, samplefreq=10,  nchains=4, savebrlens="yes", temp=0.2, path="/Applications/mrbayes-3.1.2/"){
		
	# get numbers of taxa and characters
	ntax <- length(alignment)
	nchar <- length(alignment[[1]])
	# concatenate bases
	alignment <- lapply(alignment, c2s)
	names(alignment) <- gsub(" ", "_", names(alignment))
	
	# print nexus file suitable for MrBayes
	if (mode(path)=="character") filename2 <- paste(c(path,filename),collapse="")
	filename2 <- paste(c(filename2,".bayes"), collapse="")
	write("#nexus", filename2, 		append=FALSE)
	write(" ", filename2, append=TRUE)
	write("begin data;", filename2, append=TRUE)
	write(c2s(c("\tdimensions ntax=", ntax, " nchar=", nchar, 	";")), filename2, append=TRUE)
	write("\tformat datatype=dna missing=N gap=-;", 	filename2, append=TRUE)
	write(" ", filename2, append=TRUE)
	write("matrix", filename2, append=TRUE)
	write(" ", filename2, append=TRUE)
	for (i in 1:ntax){
		s <- alignment[i]
		s <- as.character(s)
		s <- c(names(alignment[i]), toupper(s))
		write(s, filename2, append=TRUE)
		}
	write(";", filename2, append=TRUE)
	write("end;", filename2, append=TRUE)
	write(" ", filename2, append=TRUE)
	# write MrBayes block
	write("begin mrbayes;", filename2, append=TRUE)
	write(c2s(c(" lset nst=", nst, " rates=", rates, " ngammacat=", ngammacat,";")), filename2, append=TRUE)
	write(c2s(c(" mcmc nruns=", nruns, " ngen=", as.integer(ngen), " printfreq=", printfreq, " samplefreq=", samplefreq, " nchains=", nchains, " savebrlens=", savebrlens, " temp=", temp,";")), filename2, append=TRUE)
	write("end;", filename2, append=TRUE)
	
	# start mrbayes
	old.wd <- getwd()
	setwd(path)
	system(paste("./mb > execute ", filename, ".bayes", sep=""))
	setwd(old.wd)
}

