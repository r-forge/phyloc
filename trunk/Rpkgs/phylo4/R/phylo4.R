require(methods)

## base class: includes branch lengths, maybe shouldn't
setClass("phylo4",
         representation(edge="matrix",
                        edge.length="numeric",
                        Nnode="integer",
                        tip.label="character",
                        root.edge="integer"),
         prototype=list(edge=matrix(nrow=0,ncol=2),
           edge.length=numeric(0),
           Nnode=as.integer(0),
           tip.label=character(0),
           node.label=character(0),
           root.edge=as.integer(NA)),
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

## hack to ensure ape compatibility
setMethod("nTips","ANY", function(x) {
  if (class(x)=="phylo") {
    Ntip(x)
  } else stop(paste("no 'nTips' method available for",
                    deparse(substitute(x)),
                    "(class",class(x),")"))
})

setGeneric("nNodes", function(x) {
  standardGeneric("nNodes")
})
setMethod("nNodes","phylo4", function(x) {
  x@Nnode
})

setGeneric("nEdges", function(x) {
  standardGeneric("nEdges")
})
setMethod("nEdges","phylo4", function(x) {
  nrow(x@edge)
})

setGeneric("edges", function(x,order,...) {
  standardGeneric("edges")
})
setMethod("edges","phylo4", function(x,order,...) {
  x@edge
})

setGeneric("RootEdge", function(x,order,...) {
  standardGeneric("RootEdge")
})
setMethod("RootEdge","phylo4", function(x,order,...) {
  x@root.edge
})

setGeneric("isRooted", function(x) {
  standardGeneric("isRooted")
})

##  trace("isRooted",browser)
setMethod("isRooted","phylo4", function(x) {
  browser()
  !is.na(x@root.edge) ||  ## root edge explicitly defined
  ## HACK: make sure we find the right "nTips"
  tabulate(edges(x)[, 1])[phylo4::nTips(x)+1] <= 2
  ## root node (first node after last tip) has <= 2 descendants
})

setGeneric("hasEdgeLength", function(x) {
  standardGeneric("hasEdgeLength")
})
setMethod("hasEdgeLength","phylo4", function(x) {
  length(x@edge.length)>0
})

setGeneric("EdgeLength", function(x) {
  standardGeneric("EdgeLength")
})
setMethod("EdgeLength","phylo4", function(x) {
  if (!hasEdgeLength) NULL else x@edge.length
})



setGeneric("hasNodeLabels", function(x) {
  standardGeneric("hasNodeLabels")
})
setMethod("hasNodeLabels","phylo4", function(x) {
  length(x@node.label)>0
})

## labels exists already as a generic
setMethod("labels","phylo4", function(object,...) {
  object@tip.label
})

## labels exists already as a generic
setGeneric("NodeLabels", function(x) {
  standardGeneric("NodeLabels")
})
setMethod("NodeLabels","phylo4", function(object) {
  object@node.label
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
as.phylo4.phylo <- function(x,...) {
  newobj <- new("phylo4",
      edge=x$edge,
      edge.length=x$edge.length,
      Nnode=x$Nnode,
      tip.label=x$tip.label,
      root.edge=if(is.null(x$root.edge)) as.integer(NA) else x$root.edge)
  attribs = attributes(x)
  attribs$names <- NULL
  knownattr <- c("logLik","order","origin","para","xi")
  known <- names(attribs)[names(attribs) %in% knownattr]
  unknown <- names(attribs)[!names(attribs) %in% c(knownattr,"class","names")]
  if (length(unknown)>0) {
    warning(paste("unknown attributes ignored: ",unknown,collapse=" "))
  }
  for (i in known) attr(newobj,i) <- attr(x,i)
  newobj
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
setMethod("$","phylo4",function(x,name) {
  switch(name,
         edge.length=if(!hasEdgeLength(x)) NULL else x@edge.length,
         node.label=if(!hasNodeLabels(x)) NULL else x@node.label,
         root.edge=if(is.na(x@root.edge)) NULL else x@root.edge,
         attr(x,name))
})

printphylo <- function (x,printlen=6,...) {
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and", 
        nb.node, "internal nodes\n\n"))
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
    rlab <- if (isRooted(x)) "Rooted"  else "Unrooted"
    cat("\n", rlab, "; ", sep = "")
    blen <- if (hasEdgeLength(x))
        "no branch lengths"
    else "includes branch lengths"
    cat(blen, "\n", sep = "")
}


## hack for print/show 
## from http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1363.html

setGeneric("print")


setMethod("print", "phylo4", printphylo)
setMethod("show", "phylo4", function(object) printphylo(object))


## S3 generic for conversion to S4
as.phylo4 <- function (x, ...) 
{
    if (class(x) == "phylo4") 
      return(x)
    UseMethod("as.phylo4")
  }
###################################
## extensions

## extend: phylo with data
setClass("phylo4d",
         representation(nodedata="data.frame",
                        edgedata="data.frame"),
         contains="phylo4")
                       
## extend: phylo with model fit

## extend: multiPhylo4

## how does multi.phylo extend any of the other
##  (single) classes? -- can it extend more
##  than one other class?  (provided all are
##  of the same type?)
