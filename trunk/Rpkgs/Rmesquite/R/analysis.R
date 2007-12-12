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
                                           runnerMethod,
                                           returnType="D",
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
                   returnType,
                   runnerMethod,
                   module,
                   tree,
                   .jcast(categMatrix,
                          new.class = "mesquite/lib/characters/CharacterData"),
                   as.integer(charIndex));
  if (need.stop) stopMesquiteModule(module);
  result
}

mesquiteApply.TreeAndCategChar.Number <- function(mesquite=.mesquite(),
                                                  module,
                                                  tree,
                                                  categMatrix,
                                                  charIndex=1,
                                                  taxaBlock=NULL,
                                                  module.script=NULL) {
  mesquiteApply.TreeAndCategChar(mesquite,
                                 module=module,
                                 tree=tree,
                                 categMatrix=categMatrix,
                                 charIndex=charIndex,
                                 taxaBlock=taxaBlock,
                                 runnerMethod="numberForTreeAndCharacter",
                                 module.script=module.script);
}

mesquiteApply.AncestralStateCategChar <- function(mesquite=.mesquite(),
                                                  module,
                                                  tree,
                                                  categMatrix,
                                                  charIndex=1,
                                                  taxaBlock=NULL,
                                                  module.script=NULL) {
  mesquiteApply.TreeAndCategChar(mesquite,
                                 module=module,
                                 tree=tree,
                                 categMatrix=categMatrix,
                                 charIndex=charIndex,
                                 taxaBlock=taxaBlock,
                                 runnerMethod="ancestralStatesCategorical",
                                 returnType="Lmesquite/lib/characters/CharacterHistory;",
                                 module.script=module.script);
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
                                           runnerMethod="numberForTreeAndCharacter",
                                           module.script=script);
  result
}

# ==== Call's one of Mesquite's ancestral state reconstruction methods
# for categorical matrices
ancestralStatesCategorical <-  function(tree,
                                        categMatrix,
                                        charIndex=1,
                                        taxaBlock,
                                        script=NULL) {
  result <- mesquiteApply.TreeAndCategChar(module="#MargProbAncStates",
                                           tree=tree,
                                           categMatrix=categMatrix,
                                           charIndex=charIndex,
                                           taxaBlock=taxaBlock,
                                           runnerMethod="ancestralStatesCategorical",
                                           returnType="Lmesquite/lib/characters/CharacterHistory;",
                                           module.script=script);
  result
}

