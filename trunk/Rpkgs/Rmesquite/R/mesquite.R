mesquite.classpath <- function() {
}

startMesquite <- function(cp){
  if (missing(cp)) {
    cp <- mesquite.classpath();
  }
  if ((!is.null(cp)) && length(cp) > 0) {
    .jaddClassPath(cp);
  }
  mesquite.Runner <<- .jnew("mesquite/rmLink/lib/MesquiteRunner")
  .jcall(mesquite.Runner, "Lmesquite/Mesquite;", "startMesquite")
}

.mesquite <- function(mesquite) {
  if (missing(mesquite) || is.null(mesquite))
    mesquite.Runner
  else
    mesquite;
}

mesquiteReadFile <- function(mesquite,path,format=NULL){
  projectID <- .jcall(.mesquite(), "I", "readFile", path, as.character(format))
  projectID
}

getTaxaBlocks <- function(projectID){
	taxaBlocks <- .jcall(.mesquite(), "[Lmesquite/lib/Taxa;", "getTaxaBlocks", projectID)
	taxaBlocks
	
}

getNumberOfTaxa <- function(taxaBlock){
	num <- .jcall(taxaBlock, "I", "getNumTaxa")
	num
}

getTaxonNames <- function(taxaBlock){
	
}

getTreeVectors <- function(taxaBlock){
	treeVectors <- .jcall(.mesquite(), "[Lmesquite/lib/TreeVector;", "getTreeVectors", taxaBlock)
	treeVectors
}

getNumberOfTrees <- function(treeVector){
	num <- .jcall(treeVector, "I", "getNumberOfTrees")
	num
}

getPhyloTreeFromVector <- function(treeVector, treeIndex){
	treeString <- .jcall(.mesquite(), "S", "getTree", tv, as.integer(treeIndex))
	eval(parse(text = paste("tree <- ", treeString, sep="")))
	tree
}

mesquiteTaxaBlock <- function(mesquite=.mesquite(),
                              nameArray,
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

mesquiteCategoricalMatrix <- function(mesquite=.mesquite(),
                                      charMatrix,
                                      matrixName="MatrixFromR",
                                      numCols=dim(charMatrix)[2],
                                      taxaBlock,
                                      blockName=NULL) {
  if (is.character(taxaBlock)) {
    taxaBlock <- mesquiteTaxaBlock(.mesquite(), taxaBlock, blockName=blockName);
  }
  catMatrix <- .jcall(mesquite,
                      "Lmesquite/categ/lib/CategoricalData;",
                      "loadCategoricalMatrix",
                      taxaBlock,
                      matrixName,
                      as.numeric(numCols),
                      as.numeric(t(charMatrix)));
  catMatrix
}

mesquiteTree <- function(mesquite=.mesquite(),
                         newick, treeName = "TreeFromR", taxaBlock, blockName){
  if (is.character(taxaBlock)) {
    taxaBlock <- mesquiteTaxaBlock(.mesquite(), taxaBlock, blockName=blockName);
  }
  tree <- .jcall(mesquite,
                 "Lmesquite/lib/MesquiteTree;",
                 "loadTree",
                 taxaBlock,
                 treeName,
                 newick);
  tree
}

startMesquiteModule <- function(className, script = NULL){
  moduleID <- .jcall(.mesquite(),
                     "I",
                     "startModule",
                     className,
                     as.character(script));
  moduleID
}

mesquiteApply.TreeAndCharacter <- function(mesquite=.mesquite(),
                                           moduleID,
                                           tree,
                                           categMatrix,
                                           charIndex=1,
                                           taxaBlock,
                                           module.script=NULL) {
  blockName <- paste(".block.",as.integer(runif(1,min=1,max=2^31)),sep="");
  if (is.matrix(categMatrix)) {
    if (is.character(taxaBlock)) {
      taxaBlock <- mesquiteTaxaBlock(mesquite, taxaBlock, blockName);
    }
    categMatrix <- mesquiteCategoricalMatrix(mesquite,
                                             charMatrix=categMatrix,
                                             taxaBlock=taxaBlock);
  }
  if (is.character(tree)) {
    if (is.character(taxaBlock)) {
      taxaBlock <- mesquiteTaxaBlock(mesquite, taxaBlock, blockName);
    }
    tree <- mesquiteTree(mesquite, tree, taxaBlock, blockName);
  }
  if (is.character(moduleID)) {
    moduleID <- startMesquiteModule(moduleID,script=module.script);
  }
  result <- .jcall(.mesquite(),
                   "D",
                   "numberForTreeAndCharacter",
                   as.integer(moduleID),
                   mesquiteTree,
                   categMatrix,
                   as.integer(charIndex))
  result
}

biSSELikelihood <-  function(tree,
                             categMatrix,
                             charIndex=1,
                             taxaBlock,
                             script=NULL){
  result <- mesquiteApply.TreeAndCharacter(moduleID="BiSSELikelihood",
                                           tree=tree,
                                           categMatrix=categMatrix,
                                           charIndex=charIndex,
                                           taxaBlock=taxaBlock,
                                           script=script);
  result
}
