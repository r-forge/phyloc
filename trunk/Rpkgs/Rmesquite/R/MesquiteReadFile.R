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

# testing ============================================================

startMesquite(cp = c('/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/Mesquite Project/Mesquite_Folder', "/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/RMLink/Mesquite_Folder", "/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/Headless/Mesquite_Folder"));

projID <- mesquiteReadFile(path = "/Users/wmaddisn/Phylogeny & Theory/!RHackathon07/!APEmultipleTrees49LBR.nex")

taxaBlocks <- getTaxaBlocks(projectID = projID)

tb <- taxaBlocks[[1]]

treeVectors <- getTreeVectors(taxaBlock = tb)

tv <- treeVectors[[1]]

num <- getNumberOfTrees(tv)

phyloTree <- getPhyloTreeFromVector(tv, 1)

phyloTree
====================================
library(rJava)
.jinit(classpath = c('/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/Mesquite Project/Mesquite_Folder', "/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/RMLink/Mesquite_Folder", "/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/Headless/Mesquite_Folder"))
mesquiteRunner <- .jnew("mesquite/rmLink/lib/MesquiteRunner")
.jcall(mesquiteRunner, "Lmesquite/Mesquite;", "startMesquite")

projectID <- .jcall(mesquiteRunner, "I", "readFile", "/Users/wmaddisn/Phylogeny & Theory/!RHackathon07/!APEmultipleTrees49LBR.nex")
show(projectID)
taxaBlocks <- .jcall(mesquiteRunner, "[Lmesquite/lib/Taxa;", "getTaxaBlocks", projectID)
show(taxaBlocks)
tb <- taxaBlocks[[1]]
treeVectors <- .jcall(mesquiteRunner, "[Lmesquite/lib/TreeVector;", "getTreeVectors", tb)
show(treeVectors)
tv <- treeVectors[[1]]
treeString <- .jcall(mesquiteRunner, "S", "getTree", tv, as.integer(1))
eval(parse(text = paste("tree <- ", treeString, sep="")))
tree

