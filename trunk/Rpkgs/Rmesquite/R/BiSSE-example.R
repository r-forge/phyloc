startMesquite <- function(cp){
	library(rJava)
	.jinit(classpath = cp)
	mesquiteRunner <<- .jnew("mesquite/rmLink/lib/MesquiteRunner")
	.jcall(mesquiteRunner, "Lmesquite/Mesquite;", "startMesquite")
}

loadTaxaBlockToMesquite <- function(blockName = "TaxaFromR", nameArray){
	taxa <- .jcall(mesquiteRunner, "Lmesquite/lib/Taxa;", "loadTaxaBlock", blockName, nameArray)
	taxa
}

loadCategoricalMatrixToMesquite <- function(rMatrix, matrixName = "MatrixFromR", taxaBlock, numChars){
	catMatrix <- .jcall(mesquiteRunner, "Lmesquite/categ/lib/CategoricalData;", "loadCategoricalMatrix", taxaBlock, matrixName, numChars, rMatrix)
	 catMatrix
}

loadTreeToMesquite <- function(newick, treeName = "TreeFromR", taxaBlock){
	tree <- .jcall(mesquiteRunner, "Lmesquite/lib/MesquiteTree;", "loadTree", taxaBlock, treeName, newick)
	tree
}

startMesquiteModule <- function(className, script = NULL){
	moduleID <- .jcall(mesquiteRunner, "I", "startModule", className, script)
	moduleID
}

functionForTreeAndCharacter <-  function(moduleID, mesquiteTree, categMatrix, characterIndex = 1){
	result <- .jcall(mesquiteRunner, "D", "numberForTreeAndCharacter", as.integer(moduleID), mesquiteTree, categMatrix, as.integer(characterIndex))
	result
}

biSSELikelihood <-  function(mesquiteTree, categMatrix, characterIndex = 1, script = ""){
	bisse <- startMesquiteModule(className = "#BiSSELikelihood", script);
	result <- functionForTreeAndCharacter(moduleID = bisse, mesquiteTree, categMatrix, characterIndex)
	result
}

# testing ============================================================


startMesquite(cp = c('/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/Mesquite Project/Mesquite_Folder', "/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/RMLink/Mesquite_Folder", "/Users/wmaddisn/Mesquite_WorkspaceHEADLESS/Headless/Mesquite_Folder"));

taxa <- loadTaxaBlockToMesquite(nameArray = c("a1", "a2", "a3", "a4"))

matrix <- c(0,1,1,0,0,0,0,0)

catMatrix <- loadCategoricalMatrixToMesquite(rMatrix = matrix, taxaBlock=taxa, numChars=2)

tree <- loadTreeToMesquite(newick = "(((a1:1.2,a3:3.4):1.5,a4:2.8):4.2,a2):0.6;", taxaBlock = taxa)

bisse <- biSSELikelihood(mesquiteTree = tree, categMatrix = catMatrix, characterIndex = 1)
