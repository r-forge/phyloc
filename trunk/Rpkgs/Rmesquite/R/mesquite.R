#========================== Basic Mesquite functions =========================

#==== Trying to determine the location of Mesquite
mesquite.classpath <- function() {
  # default location is $HOME/Mesquite_Support_Files
  home.dir <- if (.Platform$OS.type == "windows")
    Sys.getenv("HOMEPATH") else Sys.getenv("HOME");
  mesqu.prefs.dir <- paste(home.dir,
                           .Platform$file.sep,"Mesquite_Support_Files",sep="");
  # for now the only place we have to look into is the Mesquite log,
  # so it has to have been fired up at least once, and it must have
  # been the headless version that was fired up last.
  con <- file(paste(mesqu.prefs.dir,.Platform$file.sep,"Mesquite_Log",sep=""),
              open="r");
  l <- readLines(con,n=1);
  while((length(l) > 0) && (length(grep("^Mesquite directory: ",l)) == 0)) {
    l <- readLines(con,n=1);
  }
  close(con);
  if (length(l) > 0) {
    p <- attr(regexpr("^Mesquite directory:\\s*",l),"match.length");
    return(substring(l, p+1));
  }
  character(0);
}

# ==== Starts up Mesquite.  Expensive; should be done once before
# calling Mesquite functions.
startMesquite <- function(cp){
  if (missing(cp)) {
    cp <- mesquite.classpath();
  }
  if ((!is.null(cp)) && length(cp) > 0) {
    .jaddClassPath(cp);
  }
  mesquite.Runner <<- .jnew("mesquite/rmLink/rCallsM/MesquiteRunner")
  .jcall(mesquite.Runner, "Lmesquite/Mesquite;", "startMesquite")
}

# ==== Returns the passed variable if provided and not null, and
# otherwise returns the cached instance of the Mesquite runner
.mesquite <- function(mesquite) {
  if (missing(mesquite) || is.null(mesquite))
    mesquite.Runner
  else
    mesquite;
}

# ==== Starts Mesquite module of given class.  E.g.,
# mesquite.diverse.BiSSELikelihood should be indicated by
# #BiSSELikelihood or #mesquite.diverse.BiSSELikelihood. This module
# itself may hire other modules and so forth.  The script passed can
# control the parameters of this module as well as the hiring of this
# module/s employee modules
startMesquiteModule <- function(className, script=NULL){
  if (is.null(script)) {
    script <- .jnew(class="java/lang/String");
  }
  moduleID <- .jcall(.mesquite(),
                     "Lmesquite/lib/MesquiteModule;",
                     "startModule",
                     className,
                     script);
  moduleID
}

#==== Stops Mesquite module.  Should be called when a module is no longer needed.
stopMesquiteModule <- function(module){
  .jcall(.mesquite(),
                     "V",
                     "stopModule",
                     module);
}

#==== Gets from Mesquite a list of modules available for a particular duty
availableMesquiteModules <- function(dutyName){
  moduleList <- .jcall(.mesquite(),
                     "[Lmesquite/lib/Listable;",
                     "modulesWithDuty",
                     dutyName);
  moduleList
}

#==== Shows the list of module names with explanations
showModules <- function(category){
	if (is.character(category)){
		category <- availableMesquiteModules(category)
	}
 for (i in 1:length(category))
  {
    module <- category[[i]]
    name <- .jcall(.mesquite(), "S", "getName", module)
    className <- .jcall(.mesquite(), "S", "getClassName", module)
    explanation <- .jcall(.mesquite(), "S", "getExplanation", module)
    result <- paste(name, " -- ", explanation, " (Class name: ", className, ")")
    show(result)
  }
}
#==== Shows the list of module names with explanations
showDuties <- function(){
  dutyList <- .jcall(.mesquite(),
                     "[S",
                     "dutyClasses");
  dutyList
}

#========================== Harvesting Data from Mesquite =============================
#==== Asks Mesquite to read a file.  Currently only NEXUS permitted
mesquiteReadFile <- function(mesquite,path,format="NEXUS"){
  project <- .jcall(.mesquite(), "Lmesquite/lib/MesquiteProject;", "readFile", path, as.character(format))
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


#========================== Giving Data To Mesquite =============================

#==== Creates taxa block in Mesquite with taxon names as indicated by the array
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

#==== Creates categorical matrix in Mesquite.  Requires corresponding taxa block in Mesquite to have already been created.
mesquiteCategoricalMatrix <- function(mesquite=.mesquite(),
                                      charMatrix,
                                      matrixName="MatrixFromR",
                                      numCols=dim(charMatrix)[2],
                                      taxaBlock,
                                      blockName=NULL) {
  if (is.character(taxaBlock)) {
    taxaBlock <- mesquiteTaxaBlock(mesquite,
                                   nameArray=taxaBlock, blockName=blockName);
  }
  catMatrix <- .jcall(mesquite,
                      "Lmesquite/categ/lib/CategoricalData;",
                      "loadCategoricalMatrix",
                      taxaBlock,
                      matrixName,
                      ## FIXME: change back to as.integer() on RMLink update 
                      as.integer(numCols),
                      as.integer(t(charMatrix)));
  catMatrix
}

#==== Creates tree in Mesquite from newick string. Requires corresponding taxa block in Mesquite to have already been created.
mesquiteTree <- function(mesquite=.mesquite(),
                         newick,
                         treeName="TreeFromR",
                         taxaBlock,
                         blockName=NULL){
  if (is.character(taxaBlock)) {
    taxaBlock <- mesquiteTaxaBlock(mesquite,
                                   nameArray=taxaBlock, blockName=blockName);
  }
  tree <- .jcall(mesquite,
                 "Lmesquite/lib/MesquiteTree;",
                 "loadTree",
                 taxaBlock,
                 treeName,
                 newick);
  tree
}

#========================== Calling Analyses in Mesquite =============================
#==== Calls a Mesquite module that returns values for a tree and character.
# Requires that the module be of java subclass NumberForTreeAndCharacter 
#   Does NOT that the module has already been started.  However, this means that a module is started for each request.
mesquiteApply.TreeAndCategChar <- function(mesquite=.mesquite(),
                                           module,
                                           tree,
                                           categMatrix,
                                           charIndex=1,
                                           taxaBlock=NULL,
                                           module.script=NULL) {
  blockName <- paste(".block.",as.integer(runif(1,min=1,max=2^31)),sep="");
  if (is.matrix(categMatrix)) {
    if (is.character(taxaBlock)) {
      taxaBlock <- mesquiteTaxaBlock(mesquite,
                                     nameArray=taxaBlock, blockName=blockName);
    }
    categMatrix <- mesquiteCategoricalMatrix(mesquite,
                                             charMatrix=categMatrix,
                                             taxaBlock=taxaBlock);
  }
  if (is.character(tree)) {
    if (is.character(taxaBlock)) {
      taxaBlock <- mesquiteTaxaBlock(mesquite,
                                     nameArray=taxaBlock, blockName=blockName);
    }
    tree <- mesquiteTree(mesquite, newick=tree, taxaBlock=taxaBlock);
  }
  if (is.character(module)) {
    module <- startMesquiteModule(module,script=module.script);
  }
  result <- .jcall(.mesquite(),
                   "D",
                   "numberForTreeAndCharacter",
                   module,
                   tree,
                   .jcast(categMatrix, new.class = "mesquite/lib/characters/CharacterData"),
                   as.integer(charIndex))
  #Hilmar: if the module had been newly created, need to call stopMesquiteModule
  result
}

#====  Calls Mesquite's BiSSE likelihood function.  
bisseLikelihood <-  function(tree,
                             categMatrix,
                             charIndex=1,
                             taxaBlock,
                             script=NULL){
  result <- mesquiteApply.TreeAndCategChar(module="#BiSSELikelihood",
                                           tree=tree,
                                           categMatrix=categMatrix,
                                           charIndex=charIndex,
                                           taxaBlock=taxaBlock,
                                           module.script=script);
  result
}


#==== Calls a module to reconstruct ancestral states.  Module needs to have been started.
ancestralStatesCategoricalFromModule <-  function(module, tree, categMatrix, charIndex = 1){
	#Hilmar: to be cleaned
	result <- .jcall(.mesquite(), "Lmesquite/lib/characters/CharacterHistory;", "ancestralStatesCategorical", module, tree, categMatrix, as.integer(charIndex))
	result
}

#==== Call's one of Mesquite's ancestral state reconstruction methods for categorical matrices
ancestralStatesCategorical <-  function(tree, categMatrix, charIndex = 1, script = ""){
	#Hilmar: to be cleaned
	anc <- startMesquiteModule(className = "#MargProbAncStates", script);
	result <- ancestralStatesCategoricalFromModule(module = anc, tree, categMatrix, charIndex)
	stopMesquiteModule(anc)
	result
}

