require(methods)


setClass("phylo4",
         representation(edge="matrix",
                        edge.length="numeric",
                        Nnode="integer",
                        node.label="character",
                        tip.label="character",
                        edge.label="character",
                        root.edge="integer"),
         prototype=list(edge=matrix(nrow=0,ncol=2),
           edge.length=numeric(0),
           Nnode=as.integer(0),
           tip.label=character(0),
           node.label=as.character(0),
           edge.label=as.character(0),
           ## check?
           ##           node.label = as.character(1:Nnode),
           root.edge=as.integer(NA)),
         validity=check_phylo4)

## accessor functions for all internal bits
## HORRIBLE KLUGE
nTips <- function(x,...)  { }  ## mask ape::nTips
setGeneric("nTips", function(x,...) {
  standardGeneric("nTips")
})
setMethod("nTips","phylo4", function(x,...) {
  length(x@tip.label)
})
## rm(nTips)

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
  !is.na(x@root.edge) ||  ## root edge explicitly defined
  ## HACK: make sure we find the right "nTips"
  tabulate(edges(x)[, 1])[phylo4::nTips(x)+1] <= 2
  ## root node (first node after last tip) has <= 2 descendants
  ## FIXME (?): fails with empty tree
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
  if (!hasEdgeLength(x)) NULL else x@edge.length
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
setMethod("NodeLabels","phylo4", function(x) {
  x@node.label
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
  newobj <- phylo4(x$edge, x$edge.length,
                   x$tip.label, node.label=NULL,
                   edge.label=NULL, root.edge=x$root.edge)
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

as.phylo4.multi.tree <- function(x,...) {
  newobj <- new("multiPhylo4",
                phylolist=lapply(as.phylo4.phylo,x),
                tree.names=names(x),
                tip.data=data.frame())
}

as.multi.tree.phylo4 <- function(x) {
  y <- lapply(x@phylolist,as.phylo.phylo4)
  names(y) <- x@tree.names
  if (nrow(x@tip.data)>0) warning("discarded tip data")
  class(y) <- "multi.tree"
  y
}

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
    printlen <- max(1,printlen)
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    nb.edge <- length(x$edge.label)
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and", 
        nb.node, "internal nodes\n"))

    # print tip labels
    cat("\nTip labels:\n")
    if (nb.tip > printlen) {
        cat(paste("\t", paste(x$tip.label[1:printlen], collapse = ", "), 
                  ", ...\n", sep = ""))
    } else print(x$tip.label)
    
    # print node labels
    cat("\nNode labels:\n")
    if (nb.node > printlen) {
        cat(paste("\t", paste(x$node.label[1:printlen], collapse = ", "), 
                  ", ...\n", sep = ""))
    } else print(x$node.label)
    
    # print edge labels
    cat("\nEdge labels:\n")
    if (nb.edge > printlen) {
        cat(paste("\t", paste(x$edge.label[1:printlen], collapse = ", "), 
                  ", ...\n", sep = ""))
    } else print(x$edge.label)

    # slots
    cat("\nSlots:\n")
    cat(paste("@", names(x)[1:4], sep=""),sep="\t")
    cat("\n")
    cat(paste("@", names(x)[5:7], sep=""),sep="\t")
    cat("\n")
    
    rlab <- if (isRooted(x)) "Rooted"  else "Unrooted"
    cat("\n", rlab, "; ", sep = "")
    blen <- if (hasEdgeLength(x))
        "no branch lengths"
    else "includes branch lengths"
    cat(blen, "\n\n", sep = "")
}


## hack for print/show 
## from http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1363.html

setGeneric("print")


setMethod("print", "phylo4", printphylo)
setMethod("show", "phylo4", function(object) printphylo(object))


#################
# summary phylo4
#################
## have to check that x$root.edge is NULL if missing
setMethod("summary","phylo4", function (object, quiet=FALSE)
          {
            x <- object
            res <- list()
             
            # build the result object
            res$name <- deparse(substitute(object, sys.frame(-1)))
            res$nb.tips <- length(x$tip.label)
            res$nb.nodes <- x$Nnode
              
            if(!is.null(x$edge.length)){
              res$mean.el <- mean(x$edge.length, na.rm=TRUE)
              res$var.el <- var(x$edge.length, na.rm=TRUE)
              res$sumry.el <- summary(x$edge.length)[-4]
            } else{
              res$mean.el <- NULL
              res$var.el <- NULL
              res$sumry.el <- NULL
            }
            
            res$loglik <- attr(x, "loglik")
            res$para <- attr(x, "para")
            res$xi <- attr(x, "xi")
            
            # if quiet, stop here                                        
            if(quiet) return(invisible(res))
            
            # now, print to screen is !quiet
            cat("\nPhylogenetic tree:", res$name, "\n\n")
            cat("  Number of tips:", res$nb.tips, "\n")
            cat("  Number of nodes:", res$nb.nodes, "\n")
            if(is.null(x$edge.length)) {
              cat("  No branch lengths.\n")
            } else {
              cat("  Branch lengths:\n")
              cat("    mean:", res$mean.el, "\n")
              cat("    variance:", res$var.el, "\n")
              cat("    distribution summary:\n")
              print(res$sumry.el)
            }
            
            if(!is.null(x$root.edge)){
              cat("  Root edge:", x$root.edge, "\n")
            } else {
              cat("  No root edge.\n")
            }
                       
            if (!is.null(attr(x, "loglik"))) {
              cat("Phylogeny estimated by maximum likelihood.\n")
              cat("  log-likelihood:", attr(x, "loglik"), "\n\n")
              npart <- length(attr(x, "para"))
              for (i in 1:npart) {
                cat("partition ", i, ":\n", sep = "")
                print(attr(x, "para")[[i]])
                if (i == 1)
                  next
                else cat("  contrast parameter (xi):", attr(x,"xi")[i - 1], "\n")
              }
            }
            return(invisible(res))
          } # end summary phylo4
          ) # end setMethod summary phylo4




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
         representation(tipdata="data.frame",
                        nodedata="data.frame"),
##                        edgedata="data.frame"),
         validity = function(object) {
           ## FIXME: finish this by intercepting FALSE, char string, etc.
           check_data(object)
           check_phylo4(object)
         },                   
         contains="phylo4")

setGeneric("tdata", function(x,...) {
  standardGeneric("tdata")
})
setMethod("tdata","phylo4d", function(x,which="tip",...) {
  switch(which,tip=x@tipdata,node=x@nodedata,
         allnode=rbind(x@tipdata,x@nodedata))
  ##         edge=x@edgedata)
})


phylo4d <- function(tree,data, ## edgedata=data.frame(),
                    which="tip") {
  if (which=="tip") {
    tipdata <- data
    nodedata <- data.frame()
  } else if (which=="node") {
    nodedata <- data
    tipdata <- data.frame()
  } else if (which=="all") {
    tipdata <- data[1:phylo4::nTips(tree),]
    nodedata <- data[(phylo4::nTips(tree)+1):(phylo4::nTips(tree)+nNodes(tree)),]
  }
  new("phylo4d",
      edge=tree@edge,
      edge.length=tree@edge.length,
      Nnode=tree@Nnode,
      tip.label=tree@tip.label,
      node.label=tree@node.label,
      root.edge=tree@root.edge,
      tipdata=tipdata,
      nodedata=nodedata)
##      edgedata=edgedata)
}


## I think node = internal nodes only, and edge = tips+nodes, right?
## I'm just a little confused as to why your tdata has options "tip", "allnode", and "edge"? The last two seem redundant.
## Marguerite
setMethod("summary", "phylo4d", function(object){
  x <- object
  tdata(x, "tip") -> tips
  tdata(x, "allnode") -> allnodes
  tdata(x, "edge") -> edges
  cat("Phylogenetic tree with", phylo4::nTips(x), " species and", nNodes(x), "internal nodes\n\n")
  cat("  Tree plus data object of type:", class(x), "\n")
  cat("  Species Names                :", labels(x), "\n")
  if (hasEdgeLength(x)){ 
    cat("  Has Branch Lengths (first 10):", EdgeLength(x)[1:min(length(EdgeLength(x)),10)], "\n")
  } 
  cat("  Rooted                       :", isRooted(x), "\n\n\n")
 
  cat("\nComparative data\n")
  if (nrow(tips)>0) 
    {
      cat("\nTips: data.frame with", phylo4::nTips(x), "species and", ncol(tips), "variables \n")
      print(summary(tips))
    }
  if (nrow(allnodes)>0) 
    {
      cat("\nNodes: data.frame with", nEdges(x), "species and internal nodes and", ncol(allnodes), "variables \n")                  ## May have to fix once  Node=Edge issue is settled
      print(summary(allnodes))
    }
  if (nrow(edges)>0) 
    {
      cat("\nNodes: data.frame with", nEdges(x), "internal nodes and ", ncol(edges), "variables \n")                             ## May have to fix once  Node=Edge issue is settled
      print(summary(allnodes))
    }
  
}) # end summary phylo4d


## Alternative phylo4d summary method, using phylo4 summary
## Marguerite
#setMethod("summary", "phylo4d", function(object){
#  x <- object

#  summary(as(x, "phylo4"))

#  tdata(x, "tip") -> tips
#  tdata(x, "allnode") -> allnodes
#  tdata(x, "edge") -> edges

# cat("\nComparative data\n")
# if (nrow(tips)>0) 
# {
#   cat("\nTips: data.frame with", nTips(x), "species and", ncol(tips), "variables \n")
#   print(summary(tips))
# }
# if (nrow(allnodes)>0) 
# {
#    cat("\nNodes: data.frame with", nEdges(x), "species and internal nodes and", ncol(allnodes), "variables \n")                  ## May have to fix once  Node=Edge issue is settled
# 	print(summary(allnodes))
# }
# if (nrow(edges)>0) 
# {
#    cat("\nNodes: data.frame with", nEdges(x), "internal nodes and ", ncol(edges), "variables \n")                             ## May have to fix once  Node=Edge issue is settled
# 	print(summary(allnodes))
# }

#}) # end summary phylo4d



## extend: phylo with model fit (???)
## hacked with logLik attribute from ape, but otherwise not done
  
setClass("multiPhylo4",
         representation(phylolist="list",
                        tree.names="character",
                        tipdata="data.frame"),
         contains="phylo4")


################
# show phylo4d
################
#
setMethod("show", "phylo4d", function(object){
  x <- object

  cat("\n##Comparative data##\n")
  #  print tree
  cat("\n#Tree#\n")
  printphylo(x)

  # print traits
  cat("\n#Traits#\n")
  cat("\ntipdata: data.frame containing", ncol(tdata(x,"tip")), "traits for", nrow(tdata(x,"tip")),"tips" )
  cat("\nnodedata: data.frame containing", ncol(tdata(x,"node")), "traits for", nrow(tdata(x,"node")),"nodes" )
##  cat("\nedgedata: data.frame containing", ncol(tdata(x,"edge")), "variables for", nrow(tdata(x,"edge")),"edges" )
  
}) # end summary phylo4d

## ?? setMethod("print", "phylo4", o)



################
# names methods
################
setMethod("names", signature(x = "phylo4"), function(x){
  temp <- rev(names(attributes(x)))[-1]
  return(rev(temp))
})

setMethod("names", signature(x = "phylo4d"), function(x){
  temp <- rev(names(attributes(x)))[-1]
  return(rev(temp))
})





###################
# Function .genlab
###################
# recursive function to have labels of constant length
# base = a character string
# n = number of labels
.genlab <- function(base, n) {
  f1 <- function(cha,n){
    if(nchar(cha)<n){
      cha <- paste("0",cha,sep="")
      return(f1(cha,n))
    } else {return(cha)}
  }
  w <- as.character(1:n)
  max0 <- max(nchar(w))
  w <- sapply(w, function(cha) f1(cha,max0))
  return(paste(base,w,sep=""))
}





#####################
# phylo4 constructor
#####################
#
# TEST ME . wait for validity check
#
phylo4 <- function(edge, edge.length=NULL, tip.label=NULL, node.label=NULL,
                   edge.label=NULL, root.edge=NULL,...){
  # edge
  mode(edge) <- "integer"
  if(any(is.na(edge))) stop("NA are not allowed in edge matrix")
  if(ncol(edge)>2) warning("the edge matrix has more than two columns")
  edge <- as.matrix(edge[,1:2])
  
  # edge.length
  if(!is.null(edge.length)) {
    if(!is.numeric(edge.length)) stop("edge.length is not numeric")
    edge.length <- edge.length
  } else {
    edge.length <- as.numeric(NULL)
  }

  # tip.label
  ntips <- sum(tabulate(edge[,1]) == 0)
  if(is.null(tip.label)) {
    tip.label <- .genlab("T",ntips)
  } else {
    if(length(tip.label) != ntips) stop("the tip labels are not consistent with the number of tips")
    tip.label <- as.character(tip.label)
  } 

  # node.label
  nnodes <- sum(tabulate(edge[,1]) > 0)
  if(is.null(node.label)) {
    node.label <- .genlab("N",nnodes)
  } else {
    if(length(node.label) != nnodes) stop("the node labels are not consistent with the number of nodes")
  } 

  # edge.label
  # an edge is named by the descendant
   if(is.null(edge.label)) {
     edge.label <- paste("E", edge[,2], sep="")
  } else {
    if(length(edge.label) != nrow(edge)) stop("the edge labels are not consistent with the number of edges")
     } 

  # root.edge
  if(!is.null(root.edge)) {
    if(!is.integer(root.edge)) stop("root.edge must be an integer")
    if(root.edge > nrow(edge)) stop("indicated root.edge do not exist")
  } else {
    root.edge <- as.integer(NULL)
  }
  
  # fill in the result
  res <- new("phylo4")
  res@edge <- edge
  res@edge.length <- edge.length
  res@Nnode <- nnodes
  res@tip.label <- tip.label
  res@node.label <- node.label
  res@edge.label <- edge.label
  res@root.edge <- root.edge

  if(!check_phylo4(res)) stop("Invalid object created")
  return(res)
}




######################
# phylo4d constructor
######################
#
# TEST ME . wait for validity check
phylo4d <- function(edge, edge.length=NULL, tip.label=NULL, node.label=NULL,
                   edge.label=NULL, root.edge=NULL,...){
  

}




## extend: multiPhylo4

## how does multi.phylo extend any of the other
##  (single) classes? -- can it extend more
##  than one other class?  (provided all are
##  of the same type?)
