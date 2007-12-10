### last modified CH on 13.05.2007

mrp <- function(..., file="supertree")
	{
	obj <- list(...)
	ntree <- length(obj)
    cat(paste("\nGetting matrix representation of ", ntree, " trees."))
LL <- NULL
for (k in 1:ntree)
	{
	tree <- obj[[k]]
	# get node and tip numbers
	nodes <- tree$Nnode
	tips <- length(tree$tip.label)
	# create matrix with zeros
	M <- matrix(rep(0, nodes*tips)) 
	dim(M) <- c(tips, nodes)
	
	# get MATRIX REPRESENTATION iterativly
	while (sum(M[,1])==0)
		{
		for (i in unique(tree$edge[,1])) #iteration of all internal nodes
			{
			x <- tree$edge[tree$edge[,1]==i] # node i and descending tips
			x <- x[x != i] # tips descending from node i
			if (max(x) < tips+1) 
				{
				# fill in M matrix
				for (j in 1:length(x)) M[x[j], i-tips] <- 1
				# accomodate tree$edge matrix
				next.subtending.node <- tree$edge[tree$edge[,2]==i][1]
				tree$edge[tree$edge==i] <- next.subtending.node
				elim.row <- which(tree$edge[,1] == tree$edge[,2])
				tree$edge <- tree$edge[-elim.row,]
				}
			}
		}
		
	# concatenate matices
	L <- as.list(tree$tip.label)
	names(L) <- tree$tip.label
	for (i in 1:tips)	L[[i]] <- as.character(M[i,])
	LIST <- list(LL,L)
	if (length(LL) == 0) LL <- L else LL <- concatenate.all(LIST)
	}
	
	# prepare data for nexus file
	if(is.character(file)){
		ntax <- length(LL)
		nchar <- length(LL[[1]])
		names(LL) <- gsub("-", "_", names(LL))
		filename <- paste(file, ".nex", sep="")	
		# write nexus file
		write("#nexus", filename, append=FALSE)
		write("\nbegin data;", filename, append=TRUE)
		write(paste("	dimensions ntax=", ntax, " nchar=", nchar, 	";", sep=""), filename, append=TRUE)
		write("	format datatype=standard missing=N gap=-;", 		filename, append=TRUE)
		write("\nmatrix", filename, append=TRUE)
		for (i in 1:ntax){
			s <- paste(names(LL)[i], "\t", c2s(as.character(LL[[i]])))
			write(s, filename, append=TRUE)
			}
		write(";", filename, append=TRUE)
		write("end;", filename, append=TRUE)
		cat("\nThe matrix was printed as NEXUS to your working directory.")
		}
		else invisible(LL)
}


mrp2 <- function(..., file="supertree")
	{
	obj <- list(...)
	ntree <- length(obj)
    cat(paste("\nGetting matrix representation of ", ntree, " trees."))
LL <- NULL
MAT.REP <- function(obj)
	{
	tree <- obj
	# get node and tip numbers
	nodes <- tree$Nnode
	tips <- length(tree$tip.label)
	# create matrix
	M <- matrix(rep(0, nodes*tips))
	dim(M) <- c(tips, nodes)
	
	# get MATRIX REPRESENTATION iterativly
	while (sum(M[,1])==0)
		{
		for (i in unique(tree$edge[,1]))
			{
			x <- tree$edge[tree$edge[,1]==i] # node i and descending tips
			x <- x[x != i] # tips descending from node i
			if (max(x) < tips+1) # only terminal nodes?
				{
				# fill in M matrix
				for (j in 1:length(x)) M[x[j], i-tips] <- 1
				# accomodate tree$edge matrix
				next.subtending.node <- tree$edge[tree$edge[,2]==i][1]
				tree$edge[tree$edge==i] <- next.subtending.node
				elim.row <- which(tree$edge[,1] == tree$edge[,2])
				tree$edge <- tree$edge[-elim.row,]
				}
			}
		}
		
	# concatenate matices
	L <- as.list(tree$tip.label)
	names(L) <- tree$tip.label
	for (i in 1:tips)	L[[i]] <- as.character(M[i,])
	if (length(LL) == 0) LL <- L else LL <- concatenate.all(LL, L)
	}
	LL <- lapply(obj, MAT.REP)
	
	# prepare data for nexus file
	ntax <- length(LL)
	nchar <- length(LL[[1]])
	names(LL) <- gsub("-", "_", names(LL))
	filename <- paste(file, ".nex", sep="")	
	# write nexus file
	write("#nexus", filename, append=FALSE)
	write("\nbegin data;", filename, append=TRUE)
	write(paste("	dimensions ntax=", ntax, " nchar=", nchar, 	";", sep=""), filename, append=TRUE)
	write("	format datatype=standard missing=N gap=-;", 	filename, append=TRUE)
	write("\nmatrix", filename, append=TRUE)
	for (i in 1:ntax){
		s <- paste(names(LL)[i], "\t", c2s(as.character(LL[[i]])))
		write(s, filename, append=TRUE)
		}
	write(";", filename, append=TRUE)
	write("end;", filename, append=TRUE)
	cat("\nThe matrix was printed as NEXUS to your working directory.")
}


DESCENDANTS <- function(tree, node){
	tips <- length(tree$tip.label)
	#print(paste("Looking for descendants of node ", node, sep=""))
	x <- tree$edge[,2][tree$edge[,1] == node]
	while(max(x) > tips){
		x <- x[x > tips] 
		for(h in 1:length(x)) tree$edge <- tree$edge[!tree$edge[,2] == x[h],]
		for(i in 1:length(x)) tree$edge[,1][tree$edge[,1] == x[i]] <- node
		x <- tree$edge[,2][tree$edge[,1] == node] 
		}
		invisible(x)	
	}
	
	
SISTER <- function(tree, node){
	x <- tree$edge[,1][tree$edge[,2] == node] # getmrca
	sister <- tree$edge[,2][tree$edge[,1] == x] # desc of mrca
	sister <- sister[sister != node] # eliminate node
	D <- DESCENDANTS(tree, sister) # get whole sister clade
	invisible(D)
	}
