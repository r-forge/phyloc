#========================== Calling Analyses in Mesquite =====================

# ==== Calls a Mesquite module that returns values for a tree and
# character.  Requires that the module be of java subclass
# NumberForTreeAndCharacter. Does NOT require that the module has already been
# started.  However, this means that a module is started for each
# request.
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
  need.stop <- FALSE;
  if (is.character(module)) {
    module <- startMesquiteModule(module,script=module.script);
    need.stop <- TRUE;
  }
  result <- .jcall(.mesquite(),
                   "D",
                   "numberForTreeAndCharacter",
                   module,
                   tree,
                   .jcast(categMatrix,
                          new.class = "mesquite/lib/characters/CharacterData"),
                   as.integer(charIndex));
  if (need.stop) stopMesquiteModule(module);
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

