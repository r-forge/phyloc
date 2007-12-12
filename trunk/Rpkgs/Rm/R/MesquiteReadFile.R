startMesquite <- function(cp){
	library(rJava)
	.jinit(classpath = cp)
	mesquiteRunner <<- .jnew("mesquite/rmLink/lib/MesquiteRunner")
	.jcall(mesquiteRunner, "Lmesquite/Mesquite;", "startMesquite")
}

mesquiteReadFile <- function(path){
	projectID <- .jcall(mesquiteRunner, "I", "readFile", path)
	projectID
}

getTaxaBlocks <- function(projectID){
	taxaBlocks <- .jcall(mesquiteRunner, "[Lmesquite/lib/Taxa;", "getTaxaBlocks", projectID)
	taxaBlocks
	
}
getNumberOfTaxa <- function(taxaBlock){
	num <- .jcall(taxaBlock, "I", "getNumTaxa")
	num
}
getTaxonNames <- function(taxaBlock){
	
}

getTreeVectors <- function(taxaBlock){
	treeVectors <- .jcall(mesquiteRunner, "[Lmesquite/lib/TreeVector;", "getTreeVectors", taxaBlock)
	treeVectors
}

getNumberOfTrees <- function(treeVector){
	num <- .jcall(treeVector, "I", "getNumberOfTrees")
	num
}

getPhyloTreeFromVector <- function(treeVector, treeIndex){
	treeString <- .jcall(mesquiteRunner, "S", "getTree", tv, as.integer(treeIndex))
	eval(parse(text = paste("tree <- ", treeString, sep="")))
	tree
}

