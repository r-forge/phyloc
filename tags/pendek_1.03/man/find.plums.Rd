\name{find.plums}
\alias{find.plums}
\title{ Finds all nodes in a phylogeny where all the node\sQuotes descendants are terminals }
\description{
A plum (a neologism) is a node all of whose descendants are terminals; a cherry is a plum with 
exactly two descendants.}
\usage{
find.plums(phylogeny)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{ An object of class \code{"phylo"} }
}
\details{
All cherries are plums, but not all plums are cherries.}
\value{
A list of plums, having the structure
  \item{node }{label of the plum in the phylogeny passed to this function}
  \item{members }{\code{tip.label}s of terminals in the plum}
  ...
}

\author{ Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking on the phylogeny passed to the function}



## The function is currently defined as
function(phylogeny){

#Finds the plums in a phylogeny that can contain polytomies
#A 'plum' (a neologism) is a clade composed of the descendants of a node, all of whose
#descendants are tips (i.e., a bit like cherries, but can be bigger).
#Returns a list of vectors, each vector holding the tip names from one plum

all.nodes<-unique(phylogeny$edge[,1])
nodes.with.nodes.descending<-names(table(phylogeny$edge[, 1][phylogeny$edge[, 2] < 0]))
plum.nodes<-setdiff(all.nodes,nodes.with.nodes.descending)

plum.list<-NULL
for (i in 1:length(plum.nodes))
{
	members<-phylogeny$tip.label[as.numeric(phylogeny$edge[,2][which(phylogeny$edge[,1]==plum.nodes[i])])]
	plum<-list(node=plum.nodes[i],tips=members)
	plum.list<-c(plum.list,list(plum))
  }
plum.list
  }
}
\keyword{manip}