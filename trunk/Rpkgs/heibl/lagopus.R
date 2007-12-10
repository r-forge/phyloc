#################################################################
# LAGOPUS.2.3
#					                 
# Last update 9.10.2007 by C.Heibl
# Authors: C.Heibl (core function) & N.Cusimano (MEDIAN 
# calculation + chronogram plotting).
#
# Please refer to the manual 
# available on 'www.christophheibl.de/mdtinr'
#
# Please send comments, questions, opinions, etc. 
# to 'heibl' at 'lmu.de'
#################################################################


multidivtime <- function(alignment, tree, final.tree, cons, start = "baseml", file = "LAGOPUS", noisy = 2, verbose = 1, runmode = 0, model = 3, Mgene = 0, fix_kappa = 0, kappa = 2.3, fix_alpha = 0, alpha = 1.0, Malpha = 0, ncatG = 5, fix_rho = 1, rho = 0, nparK = 0, clock = 0, nhomo = 1, getSE = 1, RateAncestor = 0, Small_Diff = 7e-6, method = 1, cleandata = 0, numgenes = 1, numsamps = 10000, sampfreq = 100, burnin = 100000, rttm = 90, brownmean = 1.11, brownsd = 1.11, minab = 1.0, newk = 0.1, othk = 0.5, thek = 0.5, bigtime = rttm*100, tipsnotcoll = 0, nodata = 0, commonbrown = 0, path.paml = FALSE, path.mdt = FALSE, screen = FALSE, transfer.files = TRUE, LogLCheck = 100)
	{
	# Check input data
	if (!is.list(alignment[[1]])) alignment <- list(alignment)
	nb.partitions <- length(alignment)
	if (!is.list(tree[[1]])){
		treelist <- vector("list", nb.partitions)
		for (i in 1:nb.partitions) treelist[[i]] <- tree
		tree <- treelist
		}
		
	cat(paste("\nFound", nb.partitions, " partition(s). Checking input data ..."))
	for (i in 1:nb.partitions){
		if (class(tree[[i]]) != "phylo") # are trees of class phylo?
        stop(paste("\nObject \"tree\" of partition", i, "is not of class \"phylo\""))
        test <-tree[[i]]$tip.label %in% names(alignment[[i]])[match(tree[[i]]$tip.label, names (alignment[[i]]))]
		if (FALSE %in% test) stop(paste("\nSome taxon names in the sequence alignment of partition", i,  "do not match the names of tips in the corresponding phylogenetic tree.\n\nPlease check your data!\n\n"))
		if (length(alignment[[i]]) != length(tree[[i]]$tip.label)) stop(paste("\nThe number of sequences in the alignment of partition", i, "does not match the numbers of tips in your phylogenetic tree.\n\nPlease check your data!\n\n"))
		}
	if (class(final.tree) != "phylo") # is final.tree of class phylo?
        stop("\nObject \"final.tree\" is not of class \"phylo\"")
    if (!is.data.frame(cons)) # is cons a data frame?
        stop("\nObject \"cons\" is not a data frame. Please refer to manual to solve this problem.")
	test <- as.vector(unlist(cons[,2:3])) %in% final.tree$tip.label
	if (!all(test)){
		misspelled <- paste(as.vector(unlist(cons[,2:3]))[!test], collapse = " ")
		stop(paste("\n\nTaxon names in the constraints do not correspond to taxon names in your final phylogenetic tree.The following taxa seem to be misspelled:\n\n", misspelled))
		}
	
	cat("\nReady to start.")
	cat(paste("\n(", date(),")", sep = ""))
	
	# fixing paths
	R.wd <- getwd()
	if (path.paml == FALSE){
		path.paml <- system("locate baseml", intern=TRUE)
		path.paml <- path.paml[1]
		path.paml <- gsub("baseml", "", path.paml)
		cat(paste("\n\npath.paml was set to", path.paml))
	}
	if (path.mdt == FALSE){
		path.mdt <- system("locate estbranches", intern=TRUE)
		path.mdt <- path.mdt[1]
		path.mdt <- gsub("estbranches", "", path.mdt)
		cat(paste("\npath.mdt was set to", path.mdt))
	}
		
	if (start == "baseml"){
		setwd(path.paml)
		for (i in 1:nb.partitions){
			ntax <- length(alignment[[i]])
			nchar <- length(alignment[[i]][[1]])
	
			## 1.step: export sequence file to PAML
	
			filename <- paste("testseq.", file, i, sep="")
			partition <- lapply(alignment[[i]], c2s)
			write(c(ntax, nchar), filename, append=FALSE)
			for (j in 1:length(partition)){
				s <- partition[j]
				s <- as.character(s)
				s <- c(names(partition[j]), s)
				write(s, filename, append=TRUE)
			}
	
			## 2.step: export tree file(s) to PAML
			filename <- paste(file, i, ".tree", sep="")
			write(paste(ntax, "1"), filename, append=FALSE)
			write.tree(tree[[i]], filename, multi.line=FALSE, append=TRUE)
	
	
			# 3.step: write baseml.ctl
			filename <- "baseml.ctl"
			write(paste("     seqfile = testseq.", file, i, sep=""), filename)
			write(paste("     outfile = paml", file, i, ".out", sep=""), filename, append=TRUE)
			write(paste("    treefile = ", file, i, ".tree\n", sep=""), filename, append=TRUE)
			write(paste("       noisy = ", noisy, "    * 0,1,2,3: how much rubbish on the screen", sep=""), filename, append=TRUE)
   			write(paste("     verbose = ", verbose, "    * 1: detailed output, 0: concise output", sep=""), filename, append=TRUE)
   			write(paste("     runmode = ", runmode, "    * 0: user tree;  1: semi-automatic;  2: automatic", sep=""), filename, append=TRUE)
    		write("                         * 3: StepwiseAddition; (4,5):PerturbationNNI", filename, append=TRUE)
    		write(paste("       model = ", model, "    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85", sep=""), filename, append=TRUE)
    		write("                         * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu", filename, append=TRUE)
    		write(paste("       Mgene = ", Mgene, "    * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff", sep=""), filename, append=TRUE)
    		write(paste("   fix_kappa = ", fix_kappa, "    * 0: estimate kappa; 1: fix kappa at value below", sep=""), filename, append=TRUE)
    		write(paste("       kappa = ", kappa, "    * initial or fixed kappa", sep=""), filename, append=TRUE)
    		write(paste("   fix_alpha = ", fix_alpha, "    * 0: estimate alpha; 1: fix alpha at value below", sep=""), filename, append=TRUE)
			write(paste("       alpha = ", alpha, "    * initial or fixed alpha, 0:infinity (constant rate)", sep=""), filename, append=TRUE)
			write(paste("      Malpha = ", Malpha, "    * 1: different alphas for genes, 0: one alpha", sep=""), filename, append=TRUE)
    		write(paste("       ncatG = ", ncatG, "    * # of categories in the dG, AdG, or nparK models of rates", sep=""), filename, append=TRUE)
    		write(paste("     fix_rho = ", fix_rho, sep=""), filename, append=TRUE)  
    		write(paste("         rho = ", rho, "    * initial or given rho,   0:no correlation", sep=""), filename, append=TRUE)
    		write(paste("       nparK = ", nparK , "    * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK", sep=""), filename, append=TRUE)
    		write(paste("       clock = ", clock, "    * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis", sep=""), filename, append=TRUE)
    		write(paste("       nhomo = ", nhomo, "    * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2", sep=""), filename, append=TRUE)
  			write(paste("       getSE = ", getSE, "    * 0: do not want them, 1: want S.E.s of estimates", sep=""), filename, append=TRUE)
 			write(paste("RateAncestor = ", RateAncestor, "    * (0,1,2): rates (alpha>0) or ancestral states", sep=""), filename, append=TRUE)
   			write(paste("  Small_Diff = ", Small_Diff, sep=""), filename, append=TRUE)
			write(paste("   cleandata = ", cleandata, "    * remove sites with ambiguity data (1:yes, 0:no)?", sep=""), filename, append=TRUE)
			write("*            ndata = 5", filename, append=TRUE)
			write("*            icode = 0    * (with RateAncestor=1. try GC in data,model=4,Mgene=4)", filename, append=TRUE)
			write("*      fix_blength = -1    * 0: ignore, -1: random, 1: initial, 2: fixed", filename, append=TRUE)
			write(paste("      method = ", method, "    * 0: simultaneous; 1: one branch at a time", sep=""), filename, append=TRUE)
			# print baseml parameters to screen
			if (screen == TRUE){
				cat(paste("\n\n'baseml.ctl' ... printed to", path.paml, sep=""))
				cat(paste("\n\n     seqfile = testseq.", file, sep=""))
				cat(paste("\n     outfile = paml", file, i, ".out", sep=""))
				cat(paste("\n    treefile = ", file, i, ".tree\n", sep=""))
				cat(paste("\n       noisy = ", noisy, "    * 0,1,2,3: how much rubbish on the screen", sep=""))
   				cat(paste("\n     verbose = ", verbose, "    * 1: detailed output, 0: concise output", sep=""))
   				cat(paste("\n     runmode = ", runmode, "    * 0: user tree;  1: semi-automatic;  2: automatic", sep=""))
    			cat("\n                    * 3: StepwiseAddition; (4,5):PerturbationNNI")
    			cat(paste("\n       model = ", model, "    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85", sep=""))
    			cat("\n                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu")
    			cat(paste("\n       Mgene = ", Mgene, "    * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff", sep=""))
    			cat(paste("\n   fix_kappa = ", fix_kappa, "    * 0: estimate kappa; 1: fix kappa at value below", sep=""))
    			cat(paste("\n       kappa = ", kappa, "  * initial or fixed kappa", sep=""))
    			cat(paste("\n   fix_alpha = ", fix_alpha, "    * 0: estimate alpha; 1: fix alpha at value below", sep=""))
				cat(paste("\n       alpha = ", alpha, "    * initial or fixed alpha, 0:infinity (constant rate)", sep=""))
				cat(paste("\n      Malpha = ", Malpha, "    * 1: different alphas for genes, 0: one alpha", sep=""))
    			cat(paste("\n       ncatG = ", ncatG, "    * # of categories in the dG, AdG, or nparK models of rates", sep=""))
    			cat(paste("\n     fix_rho = ", fix_rho, sep=""))  
    			cat(paste("\n         rho = ", rho, "    * initial or given rho,   0:no correlation", sep=""))
    			cat(paste("\n       nparK = ", nparK, "    * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK", sep=""))
    			cat(paste("\n       clock = ", clock, "    * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis", sep=""))
    			cat(paste("\n       nhomo = ", nhomo, "    * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2", sep=""))
  				cat(paste("\n       getSE = ", getSE, "    * 0: do not want them, 1: want S.E.s of estimates", sep=""))
 				cat(paste("\nRateAncestor = ", RateAncestor, "    * (0,1,2): rates (alpha>0) or ancestral states", sep=""))
   				cat(paste("\n  Small_Diff = ", Small_Diff, sep=""))
				cat(paste("\n   cleandata = ", cleandata, "    * remove sites with ambiguity data (1:yes, 0:no)?", sep=""))
				cat("\n*      ndata = 5")
				cat("\n*      icode = 0    * (with RateAncestor=1. try GC in data,model=4,Mgene=4)")
				cat("\n*fix_blength = -1   * 0: ignore, -1: random, 1: initial, 2: fixed")
				cat(paste("\n      method = ", method, "    * 0: simultaneous; 1: one branch at a time", sep=""))
			}
		
			## 5. run baseml 
			cat(paste("\n\nrunning BASEML on partition", i, "... "))
			system(paste("./baseml baseml.ctl > o", file, i, sep=""))
			cat(" finished.")
			cat(paste("\n(", date(),")", sep = ""))
	
			## 6. transfer PAML output to multidistribute directory 
			system(paste("mv paml", file, i, ".out ", path.mdt, "paml", file, i, ".out", sep=""))
		}
	}
	
	if (start == "baseml" || start == "estbranches"){
		setwd(path.mdt)	
		for (i in 1:nb.partitions){
			ntax <- length(alignment[[i]])
			nchar <- length(alignment[[i]][[1]])
	
			## 7. export sequence file to multidivtime
			filename <- "testseq"
			partition <- lapply(alignment[[i]], c2s)
			write(c(ntax, nchar), filename, append=FALSE)
			for (j in 1:length(partition)){
				s <- partition[j]
				s <- as.character(s)
				s <- c(names(partition[j]), s)
				write(s, filename, append=TRUE)
			}
		
			## 8. export tree file to multidivtime
			if (is.null(final.tree$root.edge))
				final.tree$root.edge <- 0
			ntax <- length(final.tree$tip.label)
			filename <- paste(file, ".tree", sep="")
			write(paste(ntax, "1"), filename, append=FALSE)
			final.tree$edge.length <- NULL
			write.tree(final.tree, filename, multi.line=FALSE, digits=0, append=TRUE)
	
			## 9. run paml2modelinf
			setwd(path.mdt)
			cat("\n\nrunning PAML2MODELINF ...")
			system(paste("./paml2modelinf ", "paml", file, i, ".out modelinf.", file, i, sep=""))
			cat(" finished.")
	
			## 10. write hmmcntrl.dat
			filename <- "hmmcntrl.dat"
			write("/* Which Model to use? */", filename)
			write(paste("modelinf.", file, i, sep=""), filename, append=TRUE)
			write("L  /* How much output? Options: L = Loud mode (prints more output, the", filename, append=TRUE) 
			write("      default), Q = Quiet mode (prints less output - use with parametric", filename, append=TRUE)
			write("      bootstrap) */", filename, append=TRUE)
			write("D  /* Predict Secondary Structure? Options: P= predict, D = do not predict", filename, append=TRUE) 
			write("      (the default option) */", filename, append=TRUE)
			write("N  /* Does user tree specify names (N) or specify order (O) of sequences", filename, append=TRUE) 
			write("      in sequence data file? */", filename, append=TRUE)
   			write("   /* The topology is in the file listed below */", filename, append=TRUE)
			write(paste(file,".tree", sep=""), filename, append=TRUE)   
			write("   /* End of hmmcntrl.dat */", filename, append=TRUE)
			if (screen == TRUE){
				cat(paste("\n\n'hmmcntrl.dat' ... printed to ", path.mdt, sep=""))
				cat("\n\n/* Which Model to use? */")
				cat(paste("\nmodelinf.", file, i, sep=""))
				cat("\nL  /* How much output? Options: L = Loud mode (prints more output, the") 
				cat("\n      default), Q = Quiet mode (prints less output - use with parametric")
				cat("\n      bootstrap) */")
				cat("\nD  /* Predict Secondary Structure? Options: P= predict, D = do not predict") 
				cat("\n      (the default option) */")
				cat("\nN  /* Does user tree specify names (N) or specify order (O) of sequences") 
				cat("\n      in sequence data file? */")
   				cat("\n   /* The topology is in the file listed below */")
				cat(paste("\n", file,".tree", sep=""))   
				cat("\n   /* End of hmmcntrl.dat */")
			}

			## 11. run estbranches
			cat(paste("\n\nrunning ESTBRANCHES on partition", i, "... "))
			system(paste("./estbranches oest.", file, i, " > out.oest.", file, i, sep=""))
			cat(" finished.")
			cat(paste("\n(", date(),")", sep = ""))
		}
	
		## 12. check likelihood optimization
		cat("\n\nCHECK PERFORMANCE OF LIKELIHOOD OPTIMIZATION:")
		for (i in 1:nb.partitions){
			cat(paste("\n\nPartition ", i, ":", sep=""))
			L.baseml <- scan(paste("paml", file, i, ".out", sep = ""), what="character", quiet = TRUE)
			x <- grep("lnL", L.baseml)
			L.baseml <- L.baseml[x:(x+10)]
			y <- grep("):", L.baseml)
			L.baseml <- as.numeric(L.baseml[y+1])
	
			L.estbran <- scan(paste("oest.", file, i, sep = ""), what="character", quiet = TRUE)
			z <- grep("FINAL", L.estbran) - 3 
			L.estbran <- as.numeric(L.estbran[z])
			cat(paste("\nlnL (baseml)      =", L.baseml))
			cat(paste("\nlnL (estbranches) =", L.estbran))
			p <- abs(L.baseml - L.estbran)
			if (p > LogLCheck) 
				stop(paste("\nThe difference in LogLik values exceeds threshold of ", LogLCheck, ". Exiting ...", sep=""))
			else cat(paste("\nThe difference in LogLik values does not exceed threshold of ", LogLCheck, ". Proceeding ...", sep=""))
		}
	}

	## 13. translate nodes and fix age constraints
	setwd(path.mdt)
	X <- scan(file = paste("oest.", file, "1", sep=""), what="character", quiet = TRUE)
	semikolon <- grep(";", X)
	TREE <- NULL
	for (i in 1:semikolon) TREE <- paste(TREE, X[i], sep="")
	TREE <- read.tree(text=TREE)
	# identify nodes ...
	nb.tips <- length(TREE$tip.label)
	nb.nodes <- TREE$Nnode
	nb.cons <- dim(cons)[1]
	CC <- cbind(TREE$tip.label, 1:nb.tips)
	R.nodes <- NULL
	for (i in 1:nb.cons){
		taxon1 <- as.numeric(CC[,2][CC[,1] == as.character(cons[i,2])])
		taxon2 <- as.numeric(CC[,2][CC[,1] == as.character(cons[i,3])])
		node <- getMRCA(TREE, c(taxon1, taxon2))
		R.nodes <- c(R.nodes,node)
	}
	# ... and translate them
	cat("\n\nTranslating nodes ...\n")
	mdt.nodes <- (nb.tips+nb.nodes-1):nb.tips
	r.nodes <- (nb.tips+1):(nb.tips+nb.nodes)
	node.table <- data.frame(mdt.nodes, r.nodes)
	# print node translation table
	if (screen == TRUE){
		for (i in 1:dim(node.table)[1])
			cat(paste("\nNode ", node.table[i,1], "in R corresponds to node ", node.table[i,2]), " in ESTBRANCHES/MULTIDIVTIME.", sep = "")
	}
	# Fixing age constraints
	cat("\n\nFixing age constraints ...")
	cons <- data.frame(cbind(as.character(cons$X1), R.nodes, as.character(cons$X4)))
	cons.R <- cons # store R constraints
	MDT.nodes <- NULL
	for (i in 1:nb.cons) MDT.nodes <- c(MDT.nodes, node.table[,1][node.table[,2] == cons[i,2]])
	cons <- data.frame(cbind(as.character(cons$V1), MDT.nodes, as.character(cons$V3)))
	# print control trees
	n <- length(TREE$tip.label)-1
	quartz(display = "", width = 10, height = 6)
	par(mfrow=c(1,2))
	plot(TREE, y.lim = nb.tips)
	nodelabels()
	title("Node assignment in R")
	plot(TREE, y.lim = nb.tips)
	mdt.nodes <- (nb.tips+nb.nodes-1):nb.tips
	nodelabels(mdt.nodes)
	title("Node assignment in Multidivtime")
	circle <- rep(NA, nb.nodes)
	x <- as.numeric(as.character(cons.R[,2]))-nb.tips
	circle[x] <- 21
	age <- rep(NA, nb.nodes)
	for (i in 1:nb.cons){
		y <- as.numeric(as.character(cons.R[i,2]))-nb.tips
		age[y] <- as.numeric(as.character(cons.R[i,3]))
	}
	quartz()
	plot(TREE, use.edge.length=FALSE, y.lim=nb.tips)
	nodelabels(pch=circle, cex=3.5, frame="n", col="orange", bg="orange")
	nodelabels(age, cex=0.8, frame="n")
	title("Age constraints")
	cat(" finished.")
	
	## 14. Calculating the mediam = "amount of evolution"
	cat("\n\nCalculating median of branchlength ('amount of evolution') ...")
	xx <- numeric(nb.tips + nb.nodes)
    for (i in 1:dim(TREE$edge)[1]) xx[TREE$edge[i, 2]] <- xx[TREE$edge[i,1]] + TREE$edge.length[i]
    MEDIAN<- median(xx[1:nb.tips])
    cat(" finished.")
    
	## 15. write multicntrl.dat
	filename <- paste("multicntrl.dat", sep="")
	write("/* the following lines are all needed in multicntrl.dat ...
   do not add or delete lines but change entry on left of each
   line as you see fit ...  */", filename)
	write(paste(file,".tree", sep=""), filename, append=TRUE)
	write(paste(nb.partitions, "... number of genes ... FOLLOWING LINES CONTAIN ONLY NAMES OF DATA FILES"), filename, append=TRUE)
	for (i in 1:nb.partitions)
		write(paste("oest.", file, i, sep=""), filename, append=TRUE)
	write(paste(numsamps, "... numsamps: How many times should the Markov chain be sampled?"), filename, append=TRUE)
	write(paste(sampfreq, "... sampfreq: How many cycles between samples of the Markov chain?"), filename, append=TRUE)
	write(paste(burnin, "... burnin: How many cycles before the first sample of Markov chain?"), filename, append=TRUE)
	write(paste(rttm, "... rttm: a priori expected number of time units between tip and root"), filename, append=TRUE)
	write(paste(rttm, "...  rttmsd: standard deviation of prior for time between tip and root"), filename, append=TRUE)
	write(paste(MEDIAN/rttm, "... rtrate: mean of prior distribution for rate at root node"), filename, append=TRUE)
	write(paste(MEDIAN/rttm, "... rtratesd: standard deviation of prior for rate at root node"), filename, append=TRUE)
	write(paste(brownmean, "... brownmean: mean of prior for brownian motion constant 'nu'"), filename, append=TRUE)
	write(paste(brownsd, "... brownsd: std. deviation of prior for brownian motion constant 'nu'"), filename, append=TRUE)
	write("/* the following lines are all needed (i.e., do not delete them) but you may ", filename, append=TRUE)
	write("   not want to alter entries unless you are familiar with the computer code */", filename, append=TRUE)
	write(paste(minab, "... minab: parameter for beta prior on proportional node depth"), filename, append=TRUE)
	write(paste(newk, "... newk: parameter in Markov chain proposal step"), filename, append=TRUE)
	write(paste(othk, "... othk: parameter in Markov chain proposal step"), filename, append=TRUE)
	write(paste(thek, "... thek: parameter in Markov chain proposal step"), filename, append=TRUE)
	write(paste(bigtime, "... bigtime: number higher than time units between tip and root could"), filename, append=TRUE)
	write("                   be in your wildest imagination", filename, append=TRUE)
	write("/* the program will expect the entry below to be the number of constraints", filename, append=TRUE)
	write("   and then the specified number of constraints should follow on", filename, append=TRUE)
	write("   subsequent lines */", filename, append=TRUE)
	write(paste(length(cons[,1]), "... number of constraints on node times"), filename, append=TRUE)
	for (i in 1:length(cons[,1])){
		c <- paste(cons[i,1], cons[i,2], cons[i,3])
		write(c, filename, append=TRUE)
	}
	write(paste(tipsnotcoll, "... number of tips which are not collected at time 0"), filename, append=TRUE)
	write(paste(nodata, "... nodata: 1 means approximate prior, 0 means approximate posterior"), filename, append=TRUE)
	write(paste(commonbrown, "... commonbrown: 1 if all genes have same tendency to change rate, 0 otherwise"), filename, append=TRUE)

	# 16. run multidivtime
	cat("\n\nrunning MULTIDIVTIME ...")
	system(paste("./multidivtime ", file, " > out.", file, sep=""))
	cat(" finished.")
	cat(paste("\n(", date(),")", sep = ""))

	# 17. print chronogram (Code due to N.Cusimano)
	tab<-scan(file=paste("out.", file, sep = ""), what ="char", quiet = TRUE)
	tab<-gsub(")","",tab)
	tab<-gsub(",","",tab)
	a<-grep("Actual",tab)
	t<-matrix(nrow=length(a), ncol=9)
	for(i in 3:11){
		t[,i-2]<-tab[a+i]
	}
	t<-t[, c(1,3,6,8,9)]
	mode(t)<-"numeric"
	colnames(t)<-c("node.mdt", "age", "S.D.", "min", "max")
	tree <- read.tree(paste("tree.", file, sep = ""))
	quartz(width = 14, height = 10)
	nb.tips <- length(tree$tip.label)
	nb.nodes <- tree$Nnode
	xx <- numeric(nb.tips + nb.nodes)
    for (i in 1:dim(tree$edge)[1]) xx[tree$edge[i, 2]] <- xx[tree$edge[i,1]] + tree$edge.length[i]
    m<-xx[1]+max(strwidth(tree$tip.label))
    xleft <- t[dim(t)[1],3]
	plot.phylo(tree, y.lim = nb.tips, x.lim = c(-xleft, (abs(m) + xleft*2)))
	title("Chronogram with mean and standarde deviation of node ages\nYellow bars represent 95% confidence intervals")
	axisPhylo()
	form = "line"
	L <- .last_plot.phylo$Ntip
	N <- .last_plot.phylo$Nnode
	node.ape <- c(L:1, (L+N):(L+1))
	t <- data.frame(t, node.ape)
    nodes <- (L):(L+N-1)  
    nodes2 <- t$node.ape[nodes+1]
    YY <- .last_plot.phylo$yy[nodes2]
    age <- as.numeric(as.vector(t$age[nodes+1]))
    sd <- as.numeric(as.vector(t$S.D.[nodes+1]))
    min <- as.numeric(as.vector(t$min[nodes+1]))
    max <- as.numeric(as.vector(t$max[nodes+1]))
    MAX <- max(as.numeric(as.vector(t$age)))
    min <- MAX-min
    max <- MAX-max
	x1 <- min
	x2 <- max	
	segments(x1, YY, x2, YY, lwd=12, col="yellow", cex=5)
  	nodelabels(text=paste(age,"Â±", sd, sep=" "), node=nodes2, frame="none", cex=0.9)
  	
  	# 18. transfer of files
if (transfer.files){
		setwd(R.wd)
		system(paste("mkdir", file))
		path <- paste(R.wd, "/", file, sep = "")
		setwd(path.mdt)
		for (i in 1:nb.partitions){
			system(paste("cp paml", file, i, ".out ", path, "/paml", file, i, ".out ",sep=""))
			system(paste("cp modelinf.", file, i, " " , path, "/paml", file, i, sep=""))
			system(paste("cp out.oest.", file, i, " ", path, "/out.oest.", file, i, sep=""))
			system(paste("cp oest.", file, i, " ", path, "/oest.", file, i, sep=""))
			}
		system(paste("cp ", file, ".tree ", path, "/", file, ".tree", sep=""))
		system(paste("mv tree.", file, " ", path, "/tree.", file, sep=""))
		system(paste("mv samp.", file, " ", path, "/samp.", file, sep=""))
		system(paste("mv ratio.", file, " ", path, "/ratio.", file, sep=""))
		system(paste("mv node.", file, " ", path, "/node.", file, sep=""))
		system(paste("mv out.", file, " ", path, "/out.", file, sep=""))
		cat(paste("\n\nOutput files printed to", path))
	}
	setwd(R.wd)
}

