### Last update: 31.08.2007 by CH

# Contents:
# 0. read.align
# 1. read.dna.nex
# 2. read.dna.phylip
# 3. write.dna.nex
# 4. write.dna.mrbayes
# 5. write.dna.phylip
# 6. write.dna.r8s
# 7. nex2newick
# 8. read.tree



# 0: This function reads fasta, nexus, and phylip alignment files
# Author: C.Heibl

read.align <- function(file, seq.names = NULL, tolower = FALSE)
	{
	# get complete file into memory
	X <- scan(file, what = character(), quiet 	= TRUE)
	
	
	# eliminating comments
	LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT) > 0)
    	{
    	m <- matrix(c(LEFT,RIGHT), length(LEFT), 2)
    	j <- NULL 
    	for (i in 1:length(LEFT))
    		{
    		x <- m[i,1]:m[i,2]
    		j <- c(j,x)
    		} 
    	X <- X[-j]
    	}
    	
    # get file format
    format <- "phylip"
    nexus <- grep("#nexus", X, ignore.case=TRUE)
    if (length(nexus) == 1) format <- "nexus"
    fasta <- grep(">", X, ignore.case=TRUE)
   	if (length(fasta) > 0)  format <- "fasta"
    
    if (format == "nexus")
    	{
    	# check for interleaved sequences	
    	INT <- grep("interleave", X, ignore.case=TRUE)
		interleave <- if (length(INT) > 0) TRUE	else FALSE
	
		# check for eliminate commands
		ELIM <- grep("eliminate", X, ignore.case=TRUE)
		eliminate <- if (length(ELIM) > 0) TRUE	else FALSE
	
		eliminate.begin <- grep("eliminate", X, ignore.case=TRUE)
	
  		# getting number of taxa
		ntax <- grep("ntax", X, ignore.case = TRUE) # get index ntax
		ntax <- gsub(";", "", X[ntax]) # elimintae semicolon
		ntax <- as.numeric(gsub("ntax=", "", ntax, ignore.case=TRUE))

		if (!interleave)
			{
			mat <- grep("matrix", X, ignore.case = TRUE)
			tax <- seq(mat+1, mat-1+2*ntax, by=2)
			dna <- tax+1
			obj <- as.list(X[dna])
			names(obj) <- X[tax]
			obj <- lapply(obj, s2c)
			}
		
		else
			{
			mat <- grep("matrix", X, ignore.case = TRUE)
			first.taxon <- X[mat+1]
			block.init <- grep(first.taxon, X)
			# get first block ... 
			tax <- seq(block.init[1], block.init[1]+2*ntax-2, by=2)
			dna <- tax+1
			obj <- as.list(X[dna])
			names(obj) <- X[tax]
			obj <- lapply(obj, s2c)
			# ... and add the rest
			for (i in block.init[-1])
				{
				tax <- seq(i, i+2*ntax-2, by=2)
				dna <- tax+1
				T <- as.list(X[dna])
				names(T) <- X[tax]
				T <- lapply(T, s2c)
				obj <- concatenate.all(obj, T)
				}
			}
		
		if (eliminate){
			eliminate.begin <- grep("eliminate", X, ignore.case=TRUE)
			eliminate.end <- grep("matrix", X, ignore.case=TRUE)
			eliminate <- (eliminate.begin + 1) : (eliminate.end -1)
			eliminate <- X[eliminate]
			eliminate <- gsub(";", "", eliminate)
			eliminate <- strsplit(eliminate, "-")
			eliminate <- as.integer(unlist(eliminate))
			elim.len <- length(eliminate)
			E <- vector()
			for (i in seq(1, elim.len, by = 2))
				{
				e <- eliminate[i]:eliminate[i+1]
				E <- c(E, e)
				}
			for (i in 1:ntax) obj[[i]] <- obj[[i]][-E]
			}
		}
		
	if (format == "fasta") {
		getTaxaNames <- function(x) {
        x <- sub("^ +", "", x)
        x <- sub(" +$", "", x)
        x <- sub("^['\"]", "", x)
        x <- sub("['\"]$", "", x)
        x
    	}
        start <- grep("^ {0,}>", X)
        taxa <- X[start]
        n <- length(taxa)
        obj <- vector("list", n)
        if (is.null(seq.names)) {
            taxa <- sub("^ {0,}>", "", taxa)
            seq.names <- getTaxaNames(taxa)
        }
        start <- c(start, length(X) + 1)
        for (i in 1:n) obj[[i]] <- unlist(strsplit(gsub(" ", 
            "", X[(start[i] + 1):(start[i + 1] - 1)]), NULL))
        names(obj) <- seq.names
    }

	if (format == "phylip")
    	{
    	# getting number of taxa
		ntax <- as.integer(X[1])
		tax <- seq(3, ntax*2+1, by=2)
		dna <- tax+1
		obj <- as.list(X[dna])
		names(obj) <- X[tax]
		obj <- lapply(obj, s2c)
    	}
    if (tolower) obj <- lapply(obj, tolower)
    obj
}




# 1: This function reads in a nexus file
# Author: C.Heibl

read.dna.nex <- function(filename, eliminate = TRUE)
	{
	# get complete file into memory
	X <- scan(file = filename, what = character(), quiet 	= TRUE)
	
	
	# eliminating comments
	LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT) > 0)
    	{
    	m <- matrix(c(LEFT,RIGHT), length(LEFT), 2)
    	j <- NULL 
    	for (i in 1:length(LEFT))
    		{
    		x <- m[i,1]:m[i,2]
    		j <- c(j,x)
    		} 
    	X <- X[-j]
    	}
    	
    # check for interleaved sequences	
    INT <- grep("interleave", X, ignore.case=TRUE)
	interleave <- if (length(INT) > 0) TRUE	else FALSE
	
	# check for eliminate commands
	ELIM <- grep("eliminate", X, ignore.case=TRUE)
	eliminate <- if (length(ELIM) > 0) TRUE	else FALSE
	
	eliminate.begin <- grep("eliminate", X, ignore.case=TRUE)
	
  	# getting number of taxa
	ntax <- grep("ntax", X, ignore.case = TRUE) # get index ntax
	ntax <- gsub(";", "", X[ntax]) # elimintae semicolon
	ntax <- as.numeric(gsub("ntax=", "", ntax, ignore.case=TRUE))

	
	if (!interleave)
		{
		mat <- grep("matrix", X, ignore.case = TRUE)
		tax <- seq(mat+1, mat-1+2*ntax, by=2)
		dna <- tax+1
		S <- as.list(X[dna])
		names(S) <- X[tax]
		S <- lapply(S, s2c)
		}
		
	else
		{
		mat <- grep("matrix", X, ignore.case = TRUE)
		first.taxon <- X[mat+1]
		block.init <- grep(first.taxon, X)
		# get first block ... 
		tax <- seq(block.init[1], block.init[1]+2*ntax-2, by=2)
		dna <- tax+1
		S <- as.list(X[dna])
		names(S) <- X[tax]
		S <- lapply(S, s2c)
		# ... and add the rest
		for (i in block.init[-1])
			{
			tax <- seq(i, i+2*ntax-2, by=2)
			dna <- tax+1
			T <- as.list(X[dna])
			names(T) <- X[tax]
			T <- lapply(T, s2c)
			S <- concatenate.all(S, T)
			}
		
		}
		
	if (eliminate){
		eliminate.begin <- grep("eliminate", X, ignore.case=TRUE)
		eliminate.end <- grep("matrix", X, ignore.case=TRUE)
		eliminate <- (eliminate.begin + 1) : (eliminate.end -1)
		eliminate <- X[eliminate]
		eliminate <- gsub(";", "", eliminate)
		eliminate <- strsplit(eliminate, "-")
		eliminate <- as.integer(unlist(eliminate))
		elim.len <- length(eliminate)
		E <- vector()
		for (i in seq(1, elim.len, by = 2))
			{
			e <- eliminate[i]:eliminate[i+1]
			E <- c(E, e)
			}
		for (i in 1:ntax) S[[i]] <- S[[i]][-E]
		}
	S
}


# 2: This function reads in a phylip file
# Author: C.Heibl

read.dna.phylip <- function(filename){
	# get complete file into memory
	X <- scan(file=filename, what = character(), quiet 	= TRUE)
	# getting number of taxa
	ntax <- as.integer(X[1])
	tax <- seq(3, ntax*2+1, by=2)
	dna <- tax+1
	S <- as.list(X[dna])
	names(S) <- X[tax]
	S <- lapply(S, s2c)
	#S <- gsub("?", "N", S)
}





# 3: This function writes alignments to a nexus file
# Author: C.Heibl

# datatyp and missing are set automatically
	
write.dna.nex <- function(seq, filename, interleave = FALSE)
	{
	# ADJUSTING PARAMETERS
	nb.seq <- length(seq) #get maximal sequence length
	names(seq) <- gsub(" ", "_", names(seq))
	datatype <- "standard"
	if("T" %in% seq[[1]]) datatype <- "dna"
	if("U" %in% seq[[1]]) datatype <- "rna"
	if("F" %in% seq[[1]]) datatype <- "protein"
	missing <- if("?" %in% seq[[1]]) "?" else "N"
	# FILLING UP TAXON NAMES
	y <- vector()
	for (i in 1:nb.seq) {
		name <- strsplit(names(seq[i]), "")
		y <- c(y, length(name[[1]])) 
		}
	longest.name <- max(y)
	lacking.spaces <- longest.name - y
	for (i in 1:nb.seq) {
		names(seq)[i] <- paste(names(seq)[i], c2s(rep(" ", lacking.spaces[i])))
		}
		
	# FILLING UP NUCLEOTIDES	
	x <- vector()	
	for (i in 1:nb.seq) x <- (c(x, length(seq[[i]])))
	nb.nuc <- max(x)
	#calculate lacking nucleotides
	lacking.nucleotides <- nb.nuc - x
	#create list of lacking nucleotides
	l <- as.list(lacking.nucleotides)
	# add taxon names to list of lacking nucleotides
	names(l) <- names(seq)
	# fill sequences
	for (i in 1:nb.seq) 
		{
		n <- as.numeric(l[i])
		m <- rep("-", n)
		m <- c2s(m)
		l[i] <- m
		}	
	# turn elements of l into one single string	
	l <- lapply(l, s2c)
	# append gap to original matrix
	seq.filled <- mapply(c, seq, l, SIMPLIFY=FALSE)
	
	
	# WRITE NON-INTERLEAVED NEXUS
	if (interleave == FALSE){
		
		seq.filled <- lapply(seq.filled, c2s)
	
		# write nexus file
		filename <- paste(filename, ".nex", sep="")
		write("#nexus\n", filename, append=FALSE)
		write(paste("[created on ", date(), "]\n", sep=""), 
        	filename, append = TRUE)
		write("begin data;", filename, append=TRUE)
		a <- c("	dimensions ntax=", nb.seq, " nchar=", nb.nuc, 		";")
		a <- c2s(a)
		write(a, filename, append=TRUE)
		write (paste("	format datatype=", datatype, " missing=", missing, " gap=-;", sep = ""), filename, append=TRUE)
		write("\nmatrix", filename, append=TRUE)
		for (i in 1:nb.seq){
			s <- seq.filled[[i]]
			s <- as.character(s)
			s <- paste(names(seq[i]), toupper(s), sep="  ")
			write(s, filename, append=TRUE)
		}
	}
			
	# WRITE INTERLEAVED NEXUS		
	else {
		# write nexus file
		filename <- paste(filename, ".nex", sep="")
		write("#nexus\n", filename, append=FALSE)
		write(paste("[created on ", date(), "]\n", sep=""), 
        	filename, append = TRUE)
		write("begin data;", filename, append=TRUE)
		a <- c("	dimensions ntax=", nb.seq, " nchar=", nb.nuc, 			" ;")
		a <- c2s(a)
		write(a, filename, append=TRUE)
		write (paste("	format datatype=", datatype, " missing=", missing, " gap=- interleave;", sep = ""), filename, append=TRUE)
		write("\nmatrix", filename, append=TRUE)
		
		# interleaved nucleotide matrix
		nb.blocks <- ceiling(nb.nuc/interleave)
		
		for (j in 1:nb.blocks){
			for (i in 1:nb.seq){
				s <- seq.filled[[i]]
				smin <- j*interleave-(interleave-1)
				smax <- if(j < nb.blocks) j*interleave else nb.nuc
				s <- s[smin:smax]
				s <- c2s(s)
				s <- paste(names(seq[i]), toupper(s), sep="  ")
				write(s, filename, append=TRUE)
			}
		}
	}
	write(";\n", filename, append=TRUE)
	write("end;", filename, append=TRUE)

}



# 4: This function writes DNA to a file ready-to-use with MrBayes. # Author: C.Heibl

write.dna.mrbayes <- function(alignment, filename, path="/Applications/mrbayes-3.1.2/", nst=6, rates="invgamma", ngammacat=4, nruns=2, ngen=1000000, printfreq=100, samplefreq=10,  nchains=4, savebrlens="yes", temp=0.2){
	
	# get numbers of taxa and characters
	ntax <- length(alignment)
	nchar <- length(alignment[[1]])
	# concatenate bases
	alignment <- lapply(alignment, c2s)
	names(alignment) <- gsub(" ", "_", names(alignment))
	
	# print nexus file suitible for MrBayes
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

}





# 5: This function writes DNA-Alignments into *relaxed* phylip files. It can handle unequal sequence lengthes. 
# Author: C. Heibl

write.dna.phylip <- function(seq, file, path="/Applications/RAxML-VI-HPC-2.2.3/")
	{
	# definitions:
	names(seq) <- gsub(" ", "_", names(seq))
	ntax <- length(seq)
	if (mode(path)=="character") 
		filename <- paste(path,file, ".phylip", sep="")
	else
		filename <- paste(file, ".phylip", sep="")
		
	#get maximal sequence length
	x <- vector()	
	for (i in 1: (ntax)) x <- (c(x, length(seq[[i]])))
	if (length(unique(x))!=1){
		x.max <- max(x)
		if (x.max )
			#calculate lacking nucleotides
			lacking.nucleotides <- x.max - x
			#create list of lacking nucleotides
			l <- as.list(lacking.nucleotides)
			# add taxon names to list of lacking nucleotides
			names(l) <- names(seq)

			# fill sequences
			for (i in 1:ntax) {
			n <- as.numeric(l[i])
			m <- rep("-", n)
			m <- c2s(m)
			l[i] <- m
			}
		# turn elements of l into one single string
		l <- lapply(l, s2c)
		# append gap to original matrix
		seq.filled <- mapply(c, seq, l, SIMPLIFY=F)
		seq.filled <- lapply(seq.filled, c2s)
		seq <- seq.filled
		nchar <- x.max 
	}
	
	else{
		nchar <- unique(x)
		seq <- lapply(seq, c2s)
	}
	
	# write phylip file
	write(paste(ntax, nchar), filename, append=FALSE)
	for (i in 1:length(seq)){
		s <- seq[i]
		s <- as.character(s)
		s <- c(names(seq[i]), toupper(s))
		#s <- as.factor(s)
		write(s, filename, append=TRUE)
	}
}


# !NOT READY!

write.dna.phylip2 <- function(seq, file, path="/Applications/RAxML-VI-HPC-2.2.3/")
	{
	# definitions:
	names(seq) <- gsub(" ", "_", names(seq))
	ntax <- length(seq)
	nchar <- length(seq[[1]])
	if (mode(path)=="character") 
		filename <- paste(path,file, ".phylip", sep="")
	else
		filename <- paste(file, ".phylip", sep="")

	# write phylip file
	write(paste(ntax, nchar), filename, append=FALSE)
	foo <- function(seq)
		{
		S <- paste(names(seq), as.character(c2s(seq)), sep=" ")
		write(S, filename, append=TRUE)
		}
	lapply(seq, foo)
	}




# 6: This function writes DNA to a file ready-to-use with r8s (molecular dating by penalized likelihood, M. Sanderson(2002)) 
# Author: C.Heibl

write.r8s <- function (tree, filename, path="/Applications/r8s1.71/bin/", nsites=3560, mrca=c("Cunoniaceae", "Spiraeanthemum_samoense", "Pullea_glabra", 83.5), cv=FALSE, lambda=NULL)
	{
	file <- filename
	if (mode(path)=="character") file <- paste(c		(path,file),collapse="")
	# write trees block
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""), 
        file = file, append = TRUE)
    cat("[call r8s by:\n", file = file, append = TRUE)
    cat(paste("\tcd", path, "\n"), file = file, append = TRUE)
    cat(paste("\t./r8s â€“f", filename, "-b]\n\n"), file = file, append = TRUE)
    cat("begin trees;\n", file = file, append = TRUE)
    cat("tree 1 =\n", file = file, append = TRUE)
    cat(write.tree(tree, file = "", multi.line = FALSE), 
            "\n", sep = "", file = file, append = TRUE)
    cat("end;\n\n", file = file, append = TRUE)
    # write r8s block
    cat("begin r8s;\n", file = file, append = TRUE)
    cat(paste("blformat nsites=",nsites," lengths=persite, ultrametric=no;\n"), file = file, append = TRUE)
    cat("collapse;\n", file = file, append = TRUE)
    cat(paste("prune taxon=",t$tip.label[length(t$tip.label)],";\n"), file = file, append = TRUE)
    cat(paste("mrca", mrca[1], mrca[2], mrca[3],";\n"), file = file, append = TRUE)
	cat(paste("fixage taxon=", mrca[1], "age=",mrca[4],";\n"), file = file, append = TRUE)
	if (cv==TRUE) 
		cat("divtime method=pl algorithm=tn crossv=yes cvStart=0 cvInc=0.5 cvNum=10;\n", file = file, append = TRUE)
		else
		{
		cat("divtime method=pl algorithm=tn;\n", file = file, append = TRUE)
		cat(paste("set smoothing=", lambda, ";\n"), file = file, append = TRUE)}
	cat("showage;\n", file = file, append = TRUE)
	cat("describe plot=chronogram;\n", file = file, append = TRUE)
	cat("describe plot=chrono_description;\n", file = file, append = TRUE)
	cat("bd divplot=yes;\n", file = file, append = TRUE)
    cat("end;\n", file = file, append = TRUE)
}    





write.nexus.HeiCu<-function (seqtree, file = "", translate = TRUE, 	original.data = TRUE)
	{if (class(seqtree) == "phylo")

#write Tree
	{
	cat("\n", "\t","The object you are writing into a nexus file ist a TREE", "\n", "\n")
    obj <- list(seqtree)
    if (length(obj) == 1) {
        if (class(obj[[1]]) == "phylo") 
            ntree <- 1
        else {
            obj <- unlist(obj, recursive = FALSE)
            ntree <- length(obj)
        }
    }
    else ntree <- length(obj)
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""), 
        file = file, append = TRUE)
    if (original.data) {
        if (!is.null(attr(obj[[1]], "origin"))) {
            if (!file.exists(attr(obj[[1]], "origin"))) {
                warning(paste("the file", attr(obj[[1]], "origin"), 
                  "cannot be found,\nthe original data won't be written with the tree."))
                original.data <- FALSE
            }
            else {
                ORI <- scan(file = attr(obj[[1]], "origin"), 
                  what = character(), sep = "\n", skip = 1)
                start <- grep("BEGIN TAXA;", ORI)
                ORI <- ORI[-(1:(start - 1))]
                ORI <- gsub("ENDBLOCK;", "END;", ORI)
                endblock <- grep("END;", ORI)
                start <- grep("BEGIN TREES;", ORI)
                end <- endblock[endblock > start][1]
                cat(ORI[1:(start - 1)], file = file, append = TRUE, 
                  sep = "\n")
                ORI <- ORI[-(1:end)]
            }
        }
        else original.data <- FALSE
    }
    N <- length(obj[[1]]$tip.label)
    if (!original.data) {
        cat("BEGIN TAXA;\n", file = file, append = TRUE)
        cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""), 
            file = file, append = TRUE)
        cat("\tTAXLABELS\n", file = file, append = TRUE)
        cat(paste("\t\t", obj[[1]]$tip.label, sep = ""), sep = "\n", 
            file = file, append = TRUE)
        cat("\t;\n", file = file, append = TRUE)
        cat("END;\n", file = file, append = TRUE)
    }
    cat("BEGIN TREES;\n", file = file, append = TRUE)
    if (translate) {
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        X <- paste("\t\t", 1:N, "\t", obj[[1]]$tip.label, ",", 
            sep = "")
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        token <- as.character(1:N)
        names(token) <- obj[[1]]$tip.label
        obj[[1]]$tip.label <- token
        if (ntree > 1) 
            for (i in 2:ntree) obj[[i]]$tip.label <- token[obj[[i]]$tip.label]
    }
    for (i in 1:ntree) {
        if (class(obj[[i]]) != "phylo") 
            next
        if (is.rooted(obj[[i]])) 
            cat("\tTREE * UNTITLED = [&R] ", file = file, append = TRUE)
        else cat("\tTREE * UNTITLED = [&U] ", file = file, append = TRUE)
        cat(write.tree(obj[[i]], file = "", multi.line = FALSE), 
            "\n", sep = "", file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)
    if (original.data) 
        cat(ORI, file = file, append = TRUE, sep = "\n")
	} else
	
	
#write Alignment
{
		cat("The object you are writing into a nexus file ist an alignment", "\n")
	#get maximal sequence length
	x <- vector()	
	for (i in 1: (length(seqtree))) x <- (c(x, length(seqtree[[i]])))
	x.max <- max(x)
	#calculate lacking nucleotides
	lacking.nucleotides <- x.max - x
	#create list of lacking nucleotides
	l <- as.list(lacking.nucleotides)
	# add taxon names to list of lacking nucleotides
	names(l) <- names(seqtree)

	# fill sequences
	for (i in 1:33) {
		n <- as.numeric(l[i])
		m <- rep("-", n)
		m <- c2s(m)
		l[i] <- m
		}
	# turn elements of l into one single string	
	l <- lapply(l, s2c)
	# append gap to original matrix
	seq.filled <- mapply(c, seqtree, l, SIMPLIFY=F)
	seq.filled <- lapply(seq.filled, c2s)
	# write nexus file
	write("#nexus", file, append=FALSE)
	write("begin data;", file, append=TRUE)
	a <- c("	dimensions ntax=", length(seqtree), " nchar=", x.max, 	";")
	a <- c2s(a)
	write(a, file, append=TRUE)
	write("	format datatype=dna missing=N gap=-;", 	file, append=TRUE)
	write("matrix", file, append=TRUE)
	for (i in 1:length(seqtree)){
		s <- seq.filled[i]
		s <- as.character(s)
		s <- c(names(seqtree[i]), s)
		#s <- as.factor(s)
		write(s, file, append=TRUE)
		}
	write(";", file, append=TRUE)
	write("end;", file, append=TRUE)
	}
}




write2.dna <- function (x, file, format = "interleaved", append = FALSE, nbcol = 1, colsep = " ", colw = 100, indent = NULL, blocksep = 1, nl = 40) 
{
	
    format <- match.arg(format, c("interleaved", "sequential", 
        "fasta"))
    N <- length(x)
    names(x) <- gsub(" ", "_", names(x))
    # add integers as taxon names if these are lacking
    if (is.null(names(x))) 
        names(x) <- as.character(1:N)
    if (is.null(indent)) {
        indent <- if (format %in% c("interleaved", "sequential")) 
            10
        else 0
    }
    if (indent == "") 
        indent <- 0
    if (is.numeric(indent)) 
        indent <- paste(rep(" ", indent), collapse = "")
    if (format == "interleaved") {
        if (blocksep) {
            blockseparation <- TRUE
            blocksep <- paste(rep("\n", blocksep), collapse = "")
        }
        else {
            blockseparation <- FALSE
        }
        if (nbcol < 0) 
            format <- "sequential"
    }
    zz <- if (append) 
        file(file, "a")
    else file(file, "w")
    if (format %in% c("interleaved", "sequential")) {
        S <- unique(unlist(lapply(x, length)))
        if (length(S) != 1) 
            stop("sequences must have the same length for interleaved or sequential format.")
        if (any(nchar(names(x)) > nl)) {
            warning("at least one name was longer than ", nl, " characters;\nthey will be truncated which may lead to some redundancy.\n")
            names(x) <- substr(names(x), 1, nl)
        }
        for (i in 1:N) {
            nam <- names(x)[i]
            nch <- nchar(nam)
            if (nch < 10) 
                names(x)[i] <- paste(nam, paste(rep(" ", 10 - 
                  nch), collapse = ""), sep = "")
        }
        cat(N, S, "\n", file = zz)
        if (nbcol < 0) {
            nb.block <- 1
            nbcol <- totalcol <- ceiling(S/colw)
        }
        else {
            nb.block <- ceiling(S/(colw * nbcol))
            totalcol <- ceiling(S/colw)
        }
        SEQ <- matrix(NA, N, totalcol)
        mode(SEQ) <- "character"
        for (i in 1:N) {
            X <- paste(x[[i]], collapse = "")
            for (j in 1:totalcol) SEQ[i, j] <- substr(X, 1 + 
                (j - 1) * colw, colw + (j - 1) * colw)
        }
    }
    # print phylip file
    if (format == "interleaved") {
        if (nb.block == 1) {
            for (i in 1:N) {
                cat(names(x), file = zz)
                cat(SEQ[i, ], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        }
        else {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        }
        if (nb.block > 1) {
            for (k in 2:nb.block) {
                if (blockseparation) 
                  cat(blocksep, file = zz)
                if (k == nb.block) {
                  for (i in 1:N) {
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1) * nbcol):ncol(SEQ)], 
                      sep = colsep, file = zz)
                    cat("\n", file = zz)
                  }
                }
                else {
                  for (i in 1:N) {
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1) * nbcol):(nbcol + 
                      (k - 1) * nbcol)], sep = colsep, file = zz)
                    cat("\n", file = zz)
                  }
                }
            }
        }
    }
    # print sequential phylip file
    if (format == "sequential") {
        if (nb.block == 1) {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        }
        else {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
                for (k in 2:nb.block) {
                  if (k == nb.block) {
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1) * nbcol):ncol(SEQ)], 
                      sep = colsep, file = zz)
                    cat("\n", file = zz)
                  }
                  else {
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1) * nbcol):(nbcol + 
                      (k - 1) * nbcol)], sep = colsep, file = zz)
                    cat("\n", file = zz)
                  }
                }
            }
        }
    }
    # print fasta file
    if (format == "fasta") {
        for (i in 1:N) {
            cat(">", names(x)[i], file = zz)
            cat("\n", file = zz)
            X <- paste(x[[i]], collapse = "")
            S <- length(x[[i]])
            if (nbcol < 0) {
                nb.block <- 1
                nbcol <- totalcol <- ceiling(S/colw)
            }
            else {
                totalcol <- ceiling(S/colw)
                nb.block <- ceiling(totalcol/nbcol)
            }
            SEQ <- character(totalcol)
            for (j in 1:totalcol) SEQ[j] <- substr(X, 1 + (j - 
                1) * colw, colw + (j - 1) * colw)
            for (k in 1:nb.block) {
                if (k == nb.block) {
                  cat(indent, file = zz)
                  cat(SEQ[(1 + (k - 1) * nbcol):length(SEQ)], 
                    sep = colsep, file = zz)
                  cat("\n", file = zz)
                }
                else {
                  cat(indent, file = zz)
                  cat(SEQ[(1 + (k - 1) * nbcol):(nbcol + (k - 
                    1) * nbcol)], sep = colsep, file = zz)
                  cat("\n", file = zz)
                }
            }
        }
    }
    close(zz)
}





# 7: This function transform nexus file to newickfile 
# Author: N.Cusimano
# Last update: 12.03.2007

nex2newick<-function(filename)
{
tp<-scan(file=filename, what="char")
# elimating comments
	#LEFT <- grep("\\[", tp)
    #RIGHT <- grep("\\]", tp)
    #if (length(LEFT) > 0){
   # 	m <- matrix(c(LEFT,RIGHT), length(LEFT), 2)
    #	j <- NULL 
    #	for (i in 1:length(LEFT)){
    #		M <- function(m) {x <- m[,1]:m[,2] 
    #		j <- c(j,x)
    #		return(j)}
    #		x <- m[i,1]:m[i,2]
    #		j <- c(j,x)
    #	} 
    #	tp <- tp[-j]
    #}
a1<-grep("translate", tp, ignore.case=TRUE)

if (length(a1) != 0) 
	{
	cat("\t","You want to transform a NEXUS file with translation table","\n")
#Set anchors to get translation table and trees
	a1<-a1+1
	a2<-grep("tree_", tp)
	a2<-a2[1]-2
	W<-grep("&W", tp)
	
#getting translation table
	cat("\t", "Getting translation table", "\n")
	trans<-(tp[a1:a2])
	trans<-gsub(",", "", trans)
	trans<-gsub(";", "", trans)

	no<-trans[seq(1, length(trans)-1, by=2)]
	spec<-trans[seq(2, length(trans), by=2)]
	list<-data.frame(no, spec)
	
		list$spec<-gsub("1", "ONE", list$spec)
		list$spec<-gsub("2", "TWO", list$spec)
		list$spec<-gsub("3", "THREE", list$spec)
		list$spec<-gsub("4", "FOUR", list$spec)
		list$spec<-gsub("5", "FIVE", list$spec)
		list$spec<-gsub("6", "SIX", list$spec)
		list$spec<-gsub("7", "SEVEN", list$spec)
		list$spec<-gsub("8", "EIGHT", list$spec)
		list$spec<-gsub("9", "NINE", list$spec)
		list$spec<-gsub("0", "ZERO", list$spec)
	cat("\t", "Your tree has ", length(no), " tip labels/species", "\n")

#getting trees 
	cat("\t", "You get a file with ", length(W), " trees", "\n")
	cat("\t", "Getting trees...")
	trees<-tp[W+2]
	
#substitute nodelabels with taxon names
	for (i in length(no):1)
		trees<-gsub(i, list[i,2], trees)
	
	trees<-gsub("ONE","1",trees)
	trees<-gsub("TWO","2",trees)
	trees<-gsub("THREE","3",trees)
	trees<-gsub("FOUR","4",trees)
	trees<-gsub("FIVE","5",trees)
	trees<-gsub("SIX","6",trees)
	trees<-gsub("SEVEN","7",trees)
	trees<-gsub("EIGHT","8",trees)
	trees<-gsub("NINE","9",trees)
	trees<-gsub("ZERO","0",trees)
		
		}
		
		else
	{
	cat("\t","You want to transform a NEXUS file without translation table","\n")
	W<-grep("&R", tp)
	trees<-tp[W+1]
	cat("\t", "You get a file with ", length(trees), " trees", "\n")
		}
#write Newickfile
	write(trees, file=paste(filename,".new",sep=""))
	
	cat("\t","Your NEWICK file has been printed succesfully to file", paste(filename,".new", sep=""), "\n")
	
}


# 8. This function reads in NEXUS and NEWICK trees
# Authors: Paradies and Heibl
# Last update: 14.09.07

read.phylo <- function(file = "", text = NULL, tree.names = NULL, skip = 0, comment.char = "#", ...)
	{
	# text input
	if (!is.null(text)) {
        if (!is.character(text)) stop("argument `text' must be of mode character")
        STRING <- text
        if (identical(STRING, character(0))) {
        	warning("empty character string.")
        	return(NULL)
        }
        nb.tree <- 1
        translation = FALSE
        trees <- vector("list", nb.tree)
   }
   else {
    	# file input
    	X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE)
		nexus <- grep("#nexus", X, ignore.case=TRUE)
    	nexus <- if (length(nexus) == 1) TRUE else FALSE
    	translation <- FALSE

		if (nexus){  
			LEFT <- grep("\\[", X)
    		RIGHT <- grep("\\]", X)
    		if (length(LEFT)) {
        		for (i in length(LEFT):1) {
            		if (LEFT[i] == RIGHT[i]) X[LEFT[i]] <- gsub("\\[.*\\]", "", X[LEFT[i]])
            		else {
                		X[LEFT[i]] <- gsub("\\[.*", "", X[LEFT[i]])
                		X[RIGHT[i]] <- gsub(".*\\]", "", X[RIGHT[i]])
                		if (LEFT[i] < RIGHT[i] - 1) X <- X[-((LEFT[i] + 1):(RIGHT[i] - 1))]
            		}
        		}
    		}
    		X <- gsub("ENDBLOCK;", "END;", X, ignore.case = TRUE)
    		endblock <- grep("END;", X, ignore.case = TRUE)
    		semico <- grep(";", X)
    		i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    		i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    		if (length(i2) == 1) if (i2 > i1) translation <- TRUE
    		if (translation) {
        		end <- semico[semico > i2][1]
        		x <- paste(X[i2:end], sep = "", collapse = "")
        		x <- gsub("TRANSLATE", "", x, ignore.case = TRUE)
        		x <- unlist(strsplit(x, "[,; \t]"))
        		x <- x[x != ""]
        		TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        		TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    		}
    		start <- if (translation) semico[semico > i2][1] + 1 else semico[semico > i1][1]
    		end <- endblock[endblock > i1][1] - 1
    		tree <- paste(X[start:end], sep = "", collapse = "")
    		tree <- gsub(" ", "", tree)
    		tree <- unlist(strsplit(tree, "[=;]"))
    		tree <- tree[grep("[\\(\\)]", tree)]
    		nb.tree <- length(tree)
    		STRING <- as.list(tree)
    		trees <- list()
		} ### end:nexus
		else { ### NEWICK
    		tree <- X
    		if (identical(tree, character(0))) {
        		warning("empty character string.")
        		return(NULL)
    		}
    		tree <- gsub("[ \t]", "", tree)
    		tsp <- unlist(strsplit(tree, NULL))
    		ind <- which(tsp == ";")
    		nb.tree <- length(ind) # number of trees
    		x <- c(1, ind[-nb.tree] + 1)
    		y <- ind - 1
    		if (is.na(y[1])) return(NULL)
    		else {
        		STRING <- vector("list", nb.tree)
        		for (i in 1:nb.tree) STRING[[i]] <- paste(tsp[x[i]:y[i]], 
            	sep = "", collapse = "")
    		}
    		trees <- vector("list", nb.tree)
		}
	}
	# build trees
    for (i in 1:nb.tree) {
        obj <- if (length(grep(":", STRING[[i]]))) 
            tree.build(STRING[[i]])
        else clado.build(STRING[[i]])
        if (translation) {
            for (j in 1:length(obj$tip.label)) {
                ind <- which(obj$tip.label[j] == TRANS[, 1])
                obj$tip.label[j] <- TRANS[ind, 2]
            }
            if (!is.null(obj$node.label)) {
                for (j in 1:length(obj$node.label)) {
                  ind <- which(obj$node.label[j] == TRANS[, 1])
                  obj$node.label[j] <- TRANS[ind, 2]
                }
            }
        }
        trees[[i]] <- obj
        if (sum(trees[[i]]$edge[, 1] == "-1") == 1 && dim(trees[[i]]$edge)[1] > 
            1) {
            warning("The root edge is apparently not correctly represented\nin your tree: this may be due to an extra pair of\nparentheses in your file; the returned treesect has been\ncorrected but your file may not be in a valid Newick\nformat")
            ind <- which(trees[[i]]$edge[, 1] == "-1")
            trees[[i]]$root.edge <- trees[[i]]$edge.length[ind]
            trees[[i]]$edge.length <- trees[[i]]$edge.length[-ind]
            trees[[i]]$edge <- trees[[i]]$edge[-ind, ]
            for (j in 1:length(trees[[i]]$edge)) if (as.numeric(trees[[i]]$edge[j]) < 
                0) 
                trees[[i]]$edge[j] <- as.character(as.numeric(trees[[i]]$edge[j]) + 
                  1)
            if (sum(trees[[i]]$edge[, 1] == "-1") == 1) 
                stop("There are apparently two root edges in your file: cannot read tree file")
        }
    }
    if (nb.tree == 1) 
        trees <- trees[[1]]
    else {
        names(trees) <- if (is.null(tree.names)) 
            paste("tree", 1:nb.tree, sep = "")
        else tree.names
        class(trees) <- c("multi.tree", "phylo")
    }
    if (length(grep("[\\/]", file)) == 1) 
        attr(trees, "origin") <- file
    else attr(trees, "origin") <- paste(getwd(), file, sep = "/")
    trees
} # end of read.phylo




