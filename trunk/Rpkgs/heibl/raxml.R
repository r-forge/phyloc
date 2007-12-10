# Last update: 05.02.2007

### This function gets initial categories right

raxml.start <- function (alignment, file, path="/Applications/RAxML-VI-HPC-2.2.3/")
	{
	# definitions:
	names(alignment) <- gsub(" ", "_", names(alignment))
	ntax <- length(alignment)
	nchar <- length(alignment[[1]])
	alignment <- lapply(alignment, c2s)
	filename <- paste(path, file, ".phylip", sep="")
	
	# write phylip file
	write(c(ntax, nchar), filename, append=FALSE)
	for (i in 1:length(alignment))
		{
		s <- alignment[i]
		s <- as.character(s)
		s <- c(names(alignment[i]), s)
		write(s, filename, append=TRUE)
		}
	
	# set up test calculations
	setwd(path)
	system(paste("./raxmlHPC -y -s ", file ,".phylip -m GTRCAT -n ST0", sep=""))
	system(paste("./raxmlHPC -y -s ", file ,".phylip -m GTRCAT -n ST1", sep=""))
	system(paste("./raxmlHPC -y -s ", file ,".phylip -m GTRCAT -n ST2", sep=""))
	system(paste("./raxmlHPC -y -s ", file ,".phylip -m GTRCAT -n ST3", sep=""))
	system(paste("./raxmlHPC -y -s ", file ,".phylip -m GTRCAT -n ST4", sep=""))
	system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST0 -n FI0", sep=""))
	system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST1 -n FI1", sep=""))
	system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST2 -n FI2", sep=""))
	system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST3 -n FI3", sep=""))
	system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST4 -n FI4", sep=""))
	system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST0 -n AI0", sep=""))
	system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST0 -n AI1", sep=""))
	system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST0 -n AI2", sep=""))
	system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST0 -n AI3", sep=""))
	system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phylip -t RAxML_parsimonyTree.ST0 -n AI4", sep=""))
	x <- scan("RAxML_info.FI0", what="c", sep="")
	FI0 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.FI1", what="c", sep="")
	FI1 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.FI2", what="c", sep="")
	FI2 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.FI3", what="c", sep="")
	FI3 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.FI4", what="c", sep="")
	FI4 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.AI0", what="c", sep="")
	AI0 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.AI1", what="c", sep="")
	AI1 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.AI2", what="c", sep="")
	AI2 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.AI3", what="c", sep="")
	AI3 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.AI4", what="c", sep="")
	AI4 <- x[grep("Final", x)-1]
	RESULT <- rbind(FI0, FI1, FI2, FI3, FI4, AI0, AI1, AI2, AI3, AI4)
	}
	
	# next step 1: get likelihood scores as double with six instead of three decimals
	# next step 2: clarify in precedure makes sense
	
	





### This function calls RAxML (Stamatakis 2006) from within R to estimate topology and branch lengths or do non-parametric bootstrapping. Author: C. Heibl

raxml <- function(alignment, file, runs=10, bs=FALSE, plot=FALSE, path="/Applications/RAxML-VI-HPC-2.2.3/")
	{
	# definitions:
	names(alignment) <- gsub(" ", "_", names(alignment))
	ntax <- length(alignment)
	nchar <- length(alignment[[1]])
	alignment <- lapply(alignment, c2s)
	filename <- paste(path, file, ".rax", sep="")
	
	# write phylip file
	write(c(ntax, nchar), filename, append=FALSE)
	for (i in 1:length(alignment))
		{
		s <- alignment[i]
		s <- as.character(s)
		s <- c(names(alignment[i]), s)
		write(s, filename, append=TRUE)
		}
	
	# execute search for best free	
	if (bs==FALSE)
		{
		WD <- getwd()	
		setwd(path)
		system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".rax -# ", runs, " -n ", file, ".tre", sep=""))
		x <- scan(paste(path, "RAxML_info.", file, ".tre",sep=""), what="c", sep="")
		i <-as.integer(gsub(":", "", x[grep("Best", x) + 5]))
		BEST.TREE <- read.tree(paste("RAxML_result.", file, ".tre.RUN.", i, sep=""))
		setwd(WD)
		TREENAME <- paste("tre.", file, sep="")
		write.tree(BEST.TREE, file=TREENAME)
		cat("\n\nThe best tree obtained was printed as ", TREENAME, " to your working directory", sep="")
		cat("\n\nA log file, parsimony starting trees, and the other resultant trees are stored in ", path, "\n\n",sep="")
		#if (plot==TRUE) plot(BEST.TREE)
		#invisible(BEST.TREE)
		}
		
	# do non-parametric bootstrapping	
	else
		{
		wd <- getwd()	
		setwd(path)
		system(paste("./raxmlHPC -f d -m GTRMIX -s ",file ,".rax -# ", runs," -b ", bs," -n ", file, ".boot", sep=""))
		#replicates <- read.tree(paste("RAxML_bootstrap.", file, ".BS.tre", sep=""))
		#cons <- consensus(replicates, p=0.5)
		#write.tree(cons, "con.tree")
		#system(paste("./raxmlHPC -f b -m GTRCAT -s ", file, ".phylip -z RAxML_bootstrap.", file , ".BS.tre -t con.tree -n", file, ".bipart", sep=""))
		#output.tree <- read.tree(paste("RAxML_bipartitions.", file, ".bipart", sep=""))
		#if (plot==TRUE)
			#{ 
			#plot(output.tree, file=paste(file, ".tre", sep=""))
			#nodelabels(output.tree$node.label, bg="white")
			#}
		#invisible(output.tree)	
		}
	}
	
	
	
	
	
### This functions calculates bipartions using RAxML algorithms and returns them to R. Author: C. Heibl
 
raxml.bipart <- function(seq, tree, boot, file, path="/Applications/RAxML-VI-HPC-2.2.3/")
	{
	if(!all(names(seq) %in% tree$tip.label))
		stop("Alignment and tree are not identical!")
	RWD <- getwd()
	setwd(path)
	# transfer files		
	write.dna.phylip(seq, "SEQ")
	write.tree(tree, "/Applications/RAxML-VI-HPC-2.2.3/TREE")
	write.multi.tree(boot, "/Applications/RAxML-VI-HPC-2.2.3/BOOT", multi.line=FALSE)
	
	# calculate bipartions
	system("./raxmlHPC -f b -m GTRCAT -s SEQ.phylip -z BOOT -t TREE -n BIPART")
	
	# get bipartitions
	output.tree <- read.tree("RAxML_bipartitions.BIPART")
	
	# delete files in the raxml folder
	system("rm SEQ.phylip")
	system("rm TREE")
	system("rm BOOT")
	system("rm RAxML_info.BIPART")
	system("rm RAxML_bipartitions.BIPART")
	
	setwd(RWD)
	invisible(output.tree)
	
	}
	
	
	
	
	
	
### This function writes DNA-Alignments into *relaxed* phylip files. If you want the bash commands printed to separate instructions file, set readme=TRUE and specify 'runs', 'bs.runs' and 'bs'. Author: C. Heibl

write.dna.raxml <- function(alignment, file, runs=10, bs.runs=1000, bs=14021981, read.me=TRUE, path="/Applications/RAxML-VI-HPC-2.2.3/")
	{
	# definitions:
	names(alignment) <- gsub(" ", "_", names(alignment))
	ntax <- length(alignment)
	nchar <- length(alignment[[1]])
	alignment <- lapply(alignment, c2s)
		
	# write phylip file
	write(paste(ntax, nchar), filename, append=FALSE)
	for (i in 1:length(alignment))
		{
		s <- alignment[i]
		s <- as.character(s)
		s <- c(names(alignmentalignment[i]), s)
		#s <- as.factor(s)
		write(s, filename, append=TRUE)
		}
	
	# write instructions
	if (read.me==TRUE)
		{
		filename.rm <- paste(path, "readme.", file, ".txt", sep="")
		write("# Start analysis by writing (or dragging) the following commands to your terminal shell:", filename.rm)
		write("\n\n# In order to search for the best ML tree do:", filename.rm, append=TRUE)
		write(paste("cd /", path, sep=""), filename.rm, append=TRUE)
		write(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phylip -# ", runs, " -n ", file, ".tre", sep=""), filename.rm, append=TRUE)
		write("\n\n# For inferring non-parametric bootstraps do:", filename.rm, append=TRUE)
		write(paste("cd /", path, sep=""), filename.rm, append=TRUE)
		write(paste("./raxmlHPC -f d -m GTRMIX -s ",file ,".phylip -# ", bs.runs," -b ", bs," -n ", file, ".BS.tre", sep=""), filename.rm, append=TRUE)
		}
	}





