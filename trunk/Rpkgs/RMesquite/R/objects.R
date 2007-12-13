#=================== Harvesting Data Objects from Mesquite ==================

#==== Asks Mesquite to read a file.  Currently only NEXUS permitted
mesquiteReadFile <- function(mesquite,path,format="NEXUS"){
  project <- .jcall(.mesquite(),
                    "Lmesquite/lib/MesquiteProject;",
                    "readFile",
                    path, format)
  project
}

#==== Gets taxa (OTU) blocks in project in Mesquite
getTaxaBlocks <- function(project){
	taxaBlocks <- .jcall(.mesquite(), "[Lmesquite/lib/Taxa;", "getTaxaBlocks", project)
	taxaBlocks
	
}

#==== Gets number of taxa (OTUs) in Mesquite taxa block
getNumberOfTaxa <- function(taxaBlock){
	num <- .jcall(taxaBlock, "I", "getNumTaxa")
	num
}

#==== Doesn't do anything yet
getTaxonNames <- function(taxaBlock){
	
}

#==== Gets tree vectors available for the taxa block in Mesquite
getTreeVectors <- function(taxaBlock){
	treeVectors <- .jcall(.mesquite(), "[Lmesquite/lib/TreeVector;", "getTreeVectors", taxaBlock)
	treeVectors
}

#==== Gets number of trees in Mesquite tree vector
getNumberOfTrees <- function(treeVector){
	num <- .jcall(treeVector, "I", "getNumberOfTrees")
	num
}

#==== Gets tree #treeIndex from Mesquite tree vector and returns as phylo object
getPhyloTreeFromVector <- function(treeVector, treeIndex){
	treeString <- .jcall(.mesquite(), "S", "getTree", tv, as.integer(treeIndex))
	eval(parse(text = paste("tree <- ", treeString, sep="")))
	tree
}


#==== Gets a list of Continuous matrices in Mesquite
getContinuousMatrices <- function(taxaBlock){
	#Hilmar: to be cleaned
	matrices <- .jcall(.mesquite(), "[Lmesquite/cont/lib/ContinuousData;", "getContinuousMatrices", taxaBlock)
	matrices
}

#==== Gets a list of Categorical matrices in Mesquite
getCategoricalMatrices <- function(taxaBlock){
	#Hilmar: to be cleaned
	matrices <- .jcall(.mesquite(), "[Lmesquite/categ/lib/CategoricalData;", "getCategoricalMatrices", taxaBlock)
	matrices
}

#==== Gets a list of DNA matrices in Mesquite
getDNAMatrices <- function(taxaBlock){
	#Hilmar: to be cleaned
	matrices <- .jcall(.mesquite(), "[Lmesquite/categ/lib/DNAData;", "getDNAMatrices", taxaBlock)
	matrices
}

#==== Converts Mesquite Continous matrix to R matrix
convertContinuousMatrix <- function(mesqContMatrix){
	#Hilmar: to be cleaned
	mat <- .jcall(.mesquite(), "[D", "convertMatrix", mesqContMatrix)
	rn <- getRowNames(mesqMatrix = mesqContMatrix)
	dim(mat) <- c(getNumberOfTaxa(taxaBlock = tb), getNumberOfCharacters(mesqMatrix = mesqContMatrix))
	rownames(mat) <- rn
	cn <- getColumnNames(mesqMatrix = mesqContMatrix)
	colnames(mat) <- cn
	mat
}
#==== Converts Mesquite Categorical matrix to R matrix
convertCategoricalMatrix <- function(mesqCategMatrix){
	#Hilmar: to be cleaned
	mat <- .jcall(.mesquite(), "[S", "convertMatrix", mesqCategMatrix)
	rn <- getRowNames(mesqMatrix = mesqCategMatrix)
	dim(mat) <- c(getNumberOfTaxa(taxaBlock = tb), getNumberOfCharacters(mesqMatrix = mesqCategMatrix))
	rownames(mat) <- rn
	cn <- getColumnNames(mesqMatrix = mesqCategMatrix)
	colnames(mat) <- cn
	mat
}
#==== Converts Mesquite DNA matrix to R matrix
convertDNAMatrix <- function(mesqDNAMatrix){
	#Hilmar: to be cleaned
	mat <- .jcall(.mesquite(), "[S", "convertMatrix", mesqDNAMatrix)
	rn <- getRowNames(mesqMatrix = mesqDNAMatrix)
	dim(mat) <- c(getNumberOfTaxa(taxaBlock = tb), getNumberOfCharacters(mesqMatrix = mesqDNAMatrix))
	rownames(mat) <- rn
	cn <- getColumnNames(mesqMatrix = mesqDNAMatrix)
	colnames(mat) <- cn
	mat
}

#==== Gets number of characters in Mesquite matrix
getNumberOfCharacters <- function(mesqMatrix){
	#Hilmar: to be cleaned
	num <- .jcall(mesqMatrix, "I", "getNumChars")
	num
}

#==== Gets character names from Mesquite matrix
getColumnNames <- function(mesqMatrix){
	#Hilmar: to be cleaned
	str <- .jcall(.mesquite(), "[S", "getColumnNames", .jcast(mesqMatrix, new.class = "mesquite/lib/characters/CharacterData"))
	str
}

#==== Gets taxon names from Mesquite matrix
getRowNames <- function(mesqMatrix){
	#Hilmar: to be cleaned
	str <- .jcall(.mesquite(), "[S", "getRowNames", .jcast(mesqMatrix, new.class = "mesquite/lib/characters/CharacterData"))
	str
}

#================== Converting return types from Mesquite ==================

from.RNumericMatrix <- function(obj) {
  if (class(obj) != "jobjRef") {
    stop("need to pass java object reference, not ",class(obj));
  }
  col.names <- .jcall(obj, "[Ljava/lang/String;","getColumnNames");
  row.names <- .jcall(obj, "[Ljava/lang/String;","getRowNames");
  vals <- .jcall(obj, "[[Lmesquite/lib/MesquiteNumber;","getValues");
  vals <- sapply(vals,
                 function(col) sapply(.jevalArray(col),
                                      function(obj) .jcall(obj,
                                                           "D",
                                                           "getDoubleValue")));
  if (length(row.names) < 2) {
    res <- vals;
    names(res) <- col.names;
    return(res);
  }
  matrix(vals, nrow=length(row.names), byrow=FALSE,
         dimnames=list(row.names,col.names));
}

#========================== Giving Data To Mesquite ========================

# ==== Creates taxa block in Mesquite with taxon names as indicated by
# the array
mesquiteTaxaBlock <- function(mesquite=.mesquite(),
                              nameArray=NULL,
                              blockName=NULL){
  if (is.null(blockName)) {
    blockName <- paste(".block.",as.integer(runif(1,min=1,max=2^31)),sep="");
  }
  taxa <- .jcall(mesquite,
                 "Lmesquite/lib/Taxa;",
                 "loadTaxaBlock",
                 blockName,
                 nameArray);
  taxa
}

# ==== Creates categorical matrix in Mesquite.
mesquiteCategoricalMatrix <- function(mesquite=.mesquite(),
                                      charMatrix,
                                      matrixName=NULL,
                                      numCols=dim(charMatrix)[2],
                                      taxaBlock,
                                      blockName=NULL) {
  if (is.character(taxaBlock)) {
    taxaBlock <- mesquiteTaxaBlock(mesquite,
                                   nameArray=taxaBlock, blockName=blockName);
  }
  if (is.null(matrixName)) {
    matrixName <- paste(".matrix.",as.integer(runif(1,min=1,max=2^31)),sep="");
  }
  catMatrix <- .jcall(mesquite,
                      "Lmesquite/categ/lib/CategoricalData;",
                      "loadCategoricalMatrix",
                      .jcast(taxaBlock, new.class="mesquite/lib/Taxa"),
                      matrixName,
                      as.integer(numCols),
                      as.integer(t(charMatrix)));
  catMatrix
}

# ==== Creates tree in Mesquite from newick string or from phylo object
mesquiteTree <- function(mesquite=.mesquite(),
                         tree,
                         treeName=NULL,
                         taxaBlock=NULL,
                         blockName=NULL){
  if (is.character(tree)) {
    tree <- read.tree(text=tree);
  }
  if (class(tree) != "phylo") {
    stop("tree argument must be string or class phylo, not ",class(tree),"\n");
  }
  if (is.null(taxaBlock)) {
    taxaBlock <- tree$tip.label;
  }
  if (is.character(taxaBlock)) {
    taxaBlock <- mesquiteTaxaBlock(mesquite,
                                   nameArray=taxaBlock, blockName=blockName);
  }
  if (is.null(treeName)) {
    treeName <- paste(".tree.",as.integer(runif(1,min=1,max=2^31)),sep="");
  }
  tree <- .jcall(mesquite,
                 "Lmesquite/lib/MesquiteTree;",
                 "loadTree",
                 taxaBlock,
                 treeName,
                 write.tree(tree));
  tree
}

