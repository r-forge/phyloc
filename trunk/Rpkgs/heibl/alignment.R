# Last update 26.06.2007 by NC

# Contents:
# 1. make.tiplabels
# 2. add.seqnames
# 3. add.tiplabels
# 4. concatenate.all
# 5. concatenate.matching
# 6. delete.autapo
# 7. delete.empty.cells
# 8. identical.seq

# 1: This function writes a tiplabels file to cwd.
# Suppose you got sequences by something like this:
# seq <- read.GenBank(x), then
# call make.tiplabel specifying seq and "filename".
# Author: C.Heibl

make.tiplabels <- function(seq, filename, append=FALSE){
	#fam <- rep("Oxalidaceae", length(seq))
	L <- data.frame(names(seq), attr(seq, "species"))
	write.table(L, filename, append=append, col.names=FALSE, row.names=FALSE)
	}
	

	
# 2: This function changes taxonnames.
# The index matches elements of tips$V2 according to the order of elements in tree$tip.label and account for excess of elements in in the alignment. The names of the latter are not changed.
# mode = 0 -> "Genus species", mode = 1 -> "G. species", mode = 0 -> "species"
# addinfo = c("coll", "extr", "pcr")
# Changed by CH, 26.06.07
	
add.seqnames <- function(seq, tips, mode=0, addinfo=FALSE){
	tips[[1]] <- as.character(tips[[1]]) # make character
	tips[[2]] <- as.character(tips[[2]]) # make character
	names(seq) <- gsub("_", " ", names(seq)) # eliminate "_"
	
	if(addinfo == FALSE){
	
	if (mode == 0 | mode =="Genus species"){ # "Genus species"
		for (i in 1:length(seq))
			# change only those names which are encoded in the table
			if(names(seq)[i] %in% tips$seq.id)
				names(seq)[i] <- paste(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)], "_", tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="")
		}
		
	if (mode == 1 | mode=="G. species"){ # "G. species"
		for (i in 1:length(seq))
			# change only those names which are encoded in the table
			if(names(seq)[i] %in% tips$seq.id){
				g <- strsplit(as.character(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)]), "")
				g <- unlist(g)[1]
				names(seq)[i] <- paste(g,". ", tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="")
				}
		}
		
	if (mode == 2| mode == "species"){
		for (i in 1:length(seq)){
			# change only those names which are encoded in the table
			if(names(seq)[i] %in% tips$seq.id)
				names(seq)[i] <- paste(tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)])
				
				
		}
	}
	}
	
	
	######### ADDINFO ##############
	
	else{
		if (mode == 0 | mode =="Genus species"){ # "Genus species"
		for (i in 1:length(seq)){
			# change only those names which are encoded in the table
			if(names(seq)[i] %in% tips$seq.id){
				if(addinfo == "coll"){
				coll <- paste(tips$coll[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$colln[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="")
				names(seq)[i] <- paste(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], coll, sep="_")
				}
				
				if(addinfo == "extr")
				names(seq)[i] <- paste(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$extr[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="_")
				
				if(addinfo == "pcr")
				names(seq)[i] <- paste(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$pcr[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="_")
				}
			}
		}
		
	if (mode == 1 | mode=="G. species"){ # "G. species"
		for (i in 1:length(seq))
			# change only those names which are encoded in the table
			if(names(seq)[i] %in% tips$seq.id){
				if(addinfo == "coll"){
				g <- strsplit(as.character(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)]), "")
				g <- unlist(g)[1]
				coll <- paste(tips$coll[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$colln[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="")
				names(seq)[i] <- paste(g,". ", tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], " ", coll, sep="")
				}
				
				if(addinfo == "extr"){
				g <- strsplit(as.character(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)]), "")
				g <- unlist(g)[1]
				names(seq)[i] <- paste(g,". ", tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], " ", tips$extr[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="")
				}
				
				if(addinfo == "pcr"){
				g <- strsplit(as.character(tips$gen[match(names(seq)[i], tips$seq.id, nomatch=0)]),"")
				g <- unlist(g)[1]
				names(seq)[i] <- paste(g,". ", tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], " ", tips$pcr[match(names(seq)[i], tips$seq.id, nomatch=0)],sep="")
				}
			}
		}
		
	if (mode == 2| mode == "species"){
		for (i in 1:length(seq)){
			# change only those names which are encoded in the table
			if(names(seq)[i] %in% tips$seq.id){
				
				if(addinfo == "coll"){
				coll <- paste(tips$coll[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$colln[match(names(seq)[i], tips$seq.id, nomatch=0)], sep="")
				names(seq)[i] <- paste(tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], coll, sep=" ")
				}
				
				if(addinfo=="extr")
				names(seq)[i] <- paste(tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$extr[match(names(seq)[i], tips$seq.id, nomatch=0)], sep=" ")
				
				if(addinfo=="pcr")
				names(seq)[i] <- paste(tips$spec[match(names(seq)[i], tips$seq.id, nomatch=0)], tips$pcr[match(names(seq)[i], tips$seq.id, nomatch=0)], sep=" ")
			}	
			}
		}
	}		
	return(seq)	
	}
	

	
# 3: This function changes tiplabels.
# The index matches elements of tips$V2 according to the order of elements in tree$tip.label!

add.tiplabels <- function(tree, tips){
	tree$tip.label <- tips$V2[match(tree$tip.label, tips$V1, 	nomatch=0)]
	return(tree)	# returns tree with changed tiplabels
	}

# NOTES:
# add.tiplabels can handle excess of elements in tips!,
# but not yet in tree$tip.label, try to fix nomatch!


	

# 4: This function concatenates all elements of two alignments to get a supermatrix. 
# Author: C.Heibl

concatenate.all <- function(...){
	
L <- list(...)
seq1 <- L[[1]]
for (i in 1:(length(L)-1))
	{
	seq2 <- L[[i+1]]	
	### check for multiple names
	FF1 <- factor(names(seq1))
	FF2 <- factor(names(seq2))
	if (length(levels(FF1)) != length(FF1) | length(levels(FF2)) != length(FF2)) cat("\nWARNING!\nAt least one of the alignments contains identical names, \ngenerating a problem in assigning the corresponding sequences.\n\nThe function will terminate.")

		### FIRST STEP
		# eliminate non-matching elements
		seq11 <- seq1 [names(seq1) %in% names(seq2)]
		seq22<- seq2 [names(seq2) %in% names(seq1)]
		# order elements of second alignment according to first alignment
		seq22 <- seq22[match(names(seq11), names(seq22))]
		# concatenate matching elements of two lists
		seq <- mapply(c, seq11, seq22, SIMPLIFY=FALSE)
		
		# get all sequences of seq1 that are not contained in seq2
		seq11x <- seq1 [!(names(seq1) %in% names(seq2))]
		# get all sequences of seq2 that are not contained in seq1
		seq22x <- seq2 [!(names(seq2) %in% names(seq1))]
		
		### SECOND STEP
		# add rows of "N" to taxa of seq1 that lack in seq2
		if (length(seq11x) > 0) # are there any lacking taxa in seq1
			{
			lacking.ns.11 <- c2s(rep("N", length(seq22[[1]])))
			seq111 <- lapply(seq11x, c2s)
			# create blocks of lacking nucleotides
			for (i in 1:length(seq11x)) seq111[i] <- lacking.ns.11
			seq111 <- lapply(seq111, s2c)
			seq4 <- mapply(c, seq11x, seq111, SIMPLIFY=FALSE)
			seq <- c(seq, seq4)
			}
		
		### THIRD STEP	
		# add rows of "N" for taxa that lack in seq1
		if (length(seq22x) > 0) # are there any lacking taxa in seq2
			{
			lacking.ns.22 <- c2s(rep("N", length(seq11[[1]])))
			seq222 <- lapply(seq22x, c2s)
			for (i in 1:length(seq22x)) seq222[i] <- lacking.ns.22
			seq222 <- lapply(seq222, s2c)
			seq5 <- mapply(c, seq222, seq22x, SIMPLIFY=FALSE)
			seq <- c(seq, seq5)
			seq1 <- seq
			}
		}
		return(seq)
}





# 5: This function concatenates only the matching elements of two or more alignments to get a supermatrix. 
# Author: C.Heibl
		
concatenate.matching <- function(...)
	{
	L <- list(...)
	seq1 <- L[[1]]
	for (i in 1:(length(L)-1))
		{
		seq2 <- L[[i+1]]	
		### check for multiple names
		FF1 <- factor(names(seq1))
		FF2 <- factor(names(seq2))
		if (length(levels(FF1)) != length(FF1) | length(levels(FF2)) != length(FF2)) cat("\nWARNING!\nAt least one of the alignments contains identical names, \ngenerating a problem in assigning the corresponding sequences.\n\nThe function will terminate.")

# eliminate non-matching elements
		seq1 <- seq1 [names(seq1) %in% names(seq2)]
		seq2 <- seq2 [names(seq2) %in% names(seq1)]
# order elements of second alignment according to second alignment
		seq2 <- seq2[match(names(seq1), names(seq2))]
# concatenate matching elements of two lists
		seq1 <-mapply(c, seq1, seq2, SIMPLIFY=FALSE)
		}
		invisible(seq1)	
}





### 6: This function deletes those autapomorphies from your alignment that are embedded in gaps. Bear in mind the effect on terminal branch lengths! This function should be used only in case where solely topology matters and computation is time-consuming.

delete.autapo <- function(alignment)
	{
	ntax <- length(alignment)
	nchar <- length(alignment[[1]])
	autapo <- NULL
	n <- 0
	for (i in 1:nchar)
		{
		x <- NULL
		for (j in 1:ntax) x <- c(x, alignment[[j]][i])
		y <- length(x[x=="A"]) + length(x[x=="C"]) + 
			 length(x[x=="G"]) + length(x[x=="T"]) + 
			 length(x[x=="D"]) + length(x[x=="H"]) + 
			 length(x[x=="M"]) + length(x[x=="S"]) + 
			 length(x[x=="K"]) + length(x[x=="Y"]) + 
			 length(x[x=="N"])
		if (y == 1) 
			autapo <- c(autapo, i)
			n <- n + 1	
		}
	for (i in 1:ntax) alignment[[i]] <- alignment[[i]][-autapo]
	cat(paste("\n", n, "autapomorphic sites have been deleted from the alignment."))
	return(alignment)
	}
	
	
	
	
	
# 7: This function deletes all empty colums and rows. Empty means that they contain only "-", "N", "0" and "NA". It is theresfore suited for binary and nucleotide matrices.
# Author: C.Heibl

delete.empty.cells <- function(alignment)
	{
	ntax <- length(alignment)
	nchar <- length(alignment[[1]])
	EMPTY <- NULL
	n <- 0
	
	# create vector for columns containing only N and -
	for (i in 1:nchar)
		{
		x <- NULL
		for (j in 1:ntax) x <- c(x, alignment[[j]][i])
		y <- length(x[x=="-"])  + length(x[x=="N"]) + length(x[x=="0"]) + length(x[x=="NA"])
		if (y == ntax) 
			{
			EMPTY <- c(EMPTY, i) # vector of empty columns
			n <- n + 1	
			}
		}
	# delete these columns from the alignment
	if (length(EMPTY) != 0)
		{	
		for (i in 1:ntax) alignment[[i]] <- alignment[[i]][-EMPTY]
		cat(paste("\n", n, "empty columns have been deleted from the alignment."))
		}
	else cat("\n There are no empty columns in the alignment.")
		
	EMPTY <- NULL
	m <- 0
	for (i in 1:ntax)
		{
		x <- alignment[[i]]
		y <- length(x[x=="-"])  + length(x[x=="N"]) + length(x[x=="0"]) + length(x[x=="NA"])
		if (y == nchar) 
			{
			EMPTY <- c(EMPTY, i) # vector of empty columns
			m <- m + 1	
			}
		}
		if (length(EMPTY) != 0)
			{	
			alignment <- alignment[-EMPTY]
			cat(paste("\n", m, "empty rows have been deleted from the alignment."))	
			}
		else cat("\n There are no empty rows in the alignment.")
		invisible(alignment)
	}

	
	
# 8: This function identifies identical sequences in the alignment. If delete=TRUE all but one of each set of identical sequences is deleted from the alignment
	
identical.seq <- function(alignment, delete=FALSE)
	{
	ntax <- length(alignment)
	y <- NULL
	z <- NULL
	DELETE <- NULL
	for (i in 1:(ntax-1))
		{
		for (j in (i+1):ntax)
			{
			if (identical(alignment[[i]],alignment[[j]])) 
				{
				y <- c(y, names(alignment[i]))
				z <- c(z, names(alignment[j]))
				DELETE <- unique(c(DELETE, j))
				}
			}	
		}
	DF <- data.frame(Sequence_1=y, Sequence_2=z)
	cat("\n Table of identical sequences:\n\n")
	print(DF)
	if (identical (delete, TRUE))
		{
		NAMES <- names(alignment[DELETE])
		alignment <- alignment[-DELETE]
		cat("\n ", paste(length(DELETE), " identical sequences removed from alignment:\n", sep=""))
		print(NAMES)
		}
	invisible(alignment)
	}	
