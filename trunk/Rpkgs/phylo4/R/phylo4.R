require(methods)

## base class: includes branch lengths, maybe shouldn't
setClass("phylo4",
         representation(edge="matrix",
                        edge.length="numeric",
                        Nnode="integer",
                        tip.label="character"),
         validity=function(object) {
           ## browser()
           N <- nrow(object@edge)
           if (length(object@edge.length) != N)
             return("edge lengths do not match number of edges")
           if (length(object@tip.label)+object@Nnode-1 != N)
             return("number of tip labels not consistent with number of edges and nodes")
           return(TRUE)
         })

## accessor functions for all internal bits
setGeneric("nTips", function(x) {
  standardGeneric("nTips")
})
setMethod("nTips","phylo4", function(x) {
  length(x@tip.label)
})

setGeneric("nNodes", function(x) {
  standardGeneric("nNodes")
})
setMethod("nNodes","phylo4", function(x) {
  x@Nnode
})

## labels exists already
setMethod("labels","phylo4", function(object,...) {
  object@tip.label
})

## setAs only works among S4 classes ...
## this doesn't work since phylo is not an S4 class
##   setAs("phylo","phylo4",
##         function(from) {
##           new("phylo4",
##               edge=x$edge,
##               edge.length=x$edge.length,
##               Nnode=x$Nnode,
##               tip.label=x$tip.label)
##         })
## }

## convert from phylo4 to phylo
as.phylo4.phylo <- function(x) {
  new("phylo4",
      edge=x$edge,
      edge.length=x$edge.length,
      Nnode=x$Nnode,
      tip.label=x$tip.label)
}

## convert from phylo to phylo4
as.phylo.phylo4 <- function(x) {
  y <- list(edge=x@edge,
            edge.length=x@edge.length,
            Nnode=x@Nnode,
            tip.label=x@tip.label)
  class(y) <- "phylo"
  y
}

                
## hack to allow access with $
setMethod("$","phylo4",function(x,name) attr(x,name))


printphylo <- function (x,printlen=6,...) {
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and", 
        nb.node, "internal nodes.\n\n"))
    cat("Tip labels:\n")
    if (nb.tip > printlen) {
        cat(paste("\t", paste(x$tip.label[1:printlen], collapse = ", "), 
                  ", ...\n", sep = ""))
    }
    else print(x$tip.label)
    if (!is.null(x$node.label)) {
        cat("\tNode labels:\n")
        if (nb.node > printlen) {
            cat(paste("\t", paste(x$node.label[1:printlen], collapse = ", "), 
                ",...\n", sep = ""))
        }
        else print(x$node.label)
    }
            ## HACK 2: not handling rooting properly yet.
            ##     rlab <- if (is.rooted(x)) 
            ##         "Rooted"
            ##     else "Unrooted"
            ##     cat("\n", rlab, "; ", sep = "")
    blen <- if (is.null(x$edge.length)) 
        "no branch lengths."
    else "includes branch lengths."
    cat(blen, "\n", sep = "")
}


## hack for print/show 
## from http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1363.html

setGeneric("print", function(x,...) {
  standardGeneric("print")
})


## Warning: New generic for "print" does not agree with implicit
## generic from package "base"; a new generic will be assigned with
## package ".GlobalEnv"  -- ???

setMethod("print", "phylo4", printphylo)
setMethod("show", "phylo4", function(object) printphylo(object))


## this is still hacked
setMethod("summary","phylo4",
          function(object) {
            summary(as.phylo.phylo4(object))
          })
## this doesn't QUITE work since it mangles
##  the name of the object ... should replace with summary



###################################
## extensions

## rooted/unrooted
## branch lengths/no branch lengths
## associated data/no data

## different classes, or just allow empty internal values?

## extend: phylo with data
setClass("phylod4",
         representation(data="data.frame"),
         contains="phylo4")
                       
## extend: phylo with model fit

## extend: multi.phylo

## how does multi.phylo extend any of the other
##  (single) classes? -- can it extend more
##  than one other class?  (provided all are
##  of the same type?)
