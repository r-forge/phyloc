#reads in distance matrix in Phylip format and returns dist object
#does not work with strict phylip-format labels
read.phylipdist <- function(filename){
	Dtable <- read.table(filename,row.names=1,sep="",skip=1)
	Dtable <- as.dist(Dtable)
	}


fixedLengthString <- function(s,n) {
	slength <- nchar(s)
	for(j in 1:n) s <- paste(c2s(s),"")
	s <- substr(s,1,n)
	}

write.PhylipDistMatrix <- function(distobject,dfilename,strictPhylip=FALSE){
	distmatrix <- as.matrix(distobject)
	ntaxa <- dim(distmatrix)[[1]]
	write(ntaxa,dfilename,append=FALSE)
	if (strictPhylip) {
		for(i in 1:ntaxa) {
			taxi	<- rownames(distmatrix)[[i]] 
			taxi <- fixedLengthString(taxi,10)
			rowi <- paste(taxi,distmatrix[i,1],sep="")
			for(j in 2:ntaxa) rowi <- paste(rowi,distmatrix[i,j],sep=" ")
			write(rowi,dfilename,append=TRUE)	
			}

	}
	else	
	{
	write.table(distmatrix,dfilename,append=TRUE,quote=FALSE,
				col.names=FALSE, sep="")
	}	
}

#call.Fitch <- function(dmatrix,dfilename="infile") {
#        writePhylipDistMatrix(dmatrix,dfilename)
#	writeFitchControlFile()
#	#execute Fitch
#	read.tree(treefilename)
#}


#this was to be one of several functions to write control files for
#phylip programs to be used, but Windows non-interactive 
	
write.FitchControlFile  <- function(Method="F",NoSearch=FALSE,NoUserL=FALSE,
	Power="2.0",NegBrLen=FALSE, 
	Outgroup=0,Lower=FALSE,Upper=FALSE,Subreplicates=FALSE,
	GlobalRearr=FALSE,NumberJumbles=0,NumDataSets=0,
	PrintStart=FALSE,PrintProgress=FALSE,DisplayTree=FALSE,PrintTree=TRUE,
	cfilename="fitch.ctl.txt",seed=1)
{
	write(paste("D\n",Method,sep=""),cfilename,append=FALSE)
	if (NoSearch) write("U",cfilename,append=TRUE)
	if (NoUserL) write("N",cfilename,append=TRUE)
	write(paste("P\n",Power,sep=""),cfilename,append=TRUE)
	if (NegBrLen) write("-",cfilename,append=TRUE)
	if (Outgroup) write(paste("O",Outgroup,sep=""),cfilename,append=TRUE)
	if (GlobalRearr) write("G",cfilename,append=TRUE)
	if (NumberJumbles > 0) write(paste("J\n",seed,"\n",NumberJumbles,
		sep=""),cfilename,append=TRUE)
	if (NumDataSets > 0) write(paste("M\n",NumDataSets,"\n",seed,
		sep=""),cfilename,append=TRUE)
	#have not coded for terminal output variable
	if (PrintStart) write("1",cfilename,append=TRUE)
	if (!(PrintProgress)) write("2",cfilename,append=TRUE)
	if (!(DisplayTree)) write("3",cfilename,append=TRUE)
	if (!(PrintTree)) write("4",cfilename,append=TRUE)
	write("Y",cfilename,append=TRUE)
}

deletePrefix <- function(s,n) {
	s <- s2c(s)
	for (j in 1:n) s[j] <- ""
	s <- c2s(s)
	}

read.strictphylipdist <- function(filename) {
	M <- readLines(filename)
	ntaxa <- as.integer(M[[1]])
	tlabels <- substr(M,start=1,stop=10)
	D <- matrix(NA,ntaxa,ntaxa)
	rownames(D) <- array("",ntaxa)
	for(i in 1:ntaxa) {
	rownames(D)[[i]] <- tlabels[[i+1]]
	rowi <- M[[i+1]]
	rowi <- deletePrefix(rowi,10)
	rowi <- strsplit(rowi[[1]]," ")
	rowi <- strsplit(rowi[[1]]," ") 
	for(j in 1:ntaxa) D[i,j] <- rowi[[j]]
	}
	D <- as.dist(D)
}
	
	
	
	
	

