\name{force.binary.tree}
\alias{force.binary.tree}
\title{Ensure a phylogeny contains no polytomies}
\description{
Many phylogenies contain polytomies, but many functions require bifurcating trees. This 
function replaces each polytomy by an arbitrary series of bifurcations separated by edges 
of length \code{arb.branch} (zero by default).}
\usage{
force.binary.tree(phylo, arb.branch = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{An object of class \code{"phylo"} }
  \item{arb.branch}{A number to which new branches are set. Zero is the only sensible value for 
  most purposes.}
}
\details{
If polytomies are hard, the result can be used straightforwardly in \code{pic}. If the 
polytomies are soft, use \code{pic.lm} with the polytomous phylogeny. The function \code{map.nodes} 
maps nodes between the polytomous and bifurcating trees.}
\value{
An object of class \code{"phylo"} with no polytomies.
}

\author{David Orme}
\note{Modified by Andy Purvis (to fix the internal structure of the tree) and Rich 
Grenyer (to include the multi-line=FALSE directive)}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\seealso{\code{\link{map.nodes}}}


## The function is currently defined as
function(phylo, arb.branch =0){

# Written by David Orme
# Takes a polytomous tree and arbitrarily resolves the polytomies
# into a sequence of dichotomies with edge length equal to arb.branch
# I can't see any reason for this to be anything other than zero but
# using no-zero values was helpful in debugging the function.

# Minor change by Andy Purvis
# Now writes the tree out and reads it back in to give it standard ape phylo structure
# Previously, only some ape commands recognised the binary tree as valid

# Minor modification by Rich Grenyer to include new nulti-line=F directive 

if(class(phylo) != "phylo") stop(substitute(phylo), " is not of class 'phylo'.")

if(is.binary.tree(phylo)){
	cat(substitute(phylo), " is already binary.")
	invisible(phylo)
  } else {	
	# get vector of polytomous node ids
	polytomies <- names(which(table(phylo$edge[,1]) > 2))
	
	# initial next node number
	new.node <-  min(as.numeric(phylo$edge)) -1
		
	# loop through the polytomies:
	for(current.polytomy in polytomies){
	
		# store the original current polytomy and
		# its parent for later repair
		orig.polytomy <- current.polytomy
		orig.parent <- phylo$edge[,1][which(phylo$edge[,2] == orig.polytomy)]
		
		# get the indices of the edges forming the polytomy
		poly.edges <- which(phylo$edge[,1] == current.polytomy)
		
		# loop over the edges after the first two
		# use indices so that poly.edges can
		# be incremented within the loop
		
		for(vals in 3:length(poly.edges)){
			# get the current edge
			current.edge <-  poly.edges[vals]
			# insert a new line in edges above the current edge
			# containing the values:
			# new.node   current.polytomy
			phylo$edge <- rbind(phylo$edge[1:(current.edge-1),], # top half
				c(new.node, current.polytomy) # new vals
				,phylo$edge[current.edge:dim(phylo$edge)[1],] )# bottom half
						
			# connect up the newly created bifurcation
			phylo$edge[current.edge+1,1] <- new.node
			
			# insert the new length edge into edge.length
			phylo$edge.length <- c(phylo$edge.length[1:current.edge -1],
				arb.branch, 
phylo$edge.length[current.edge:length(phylo$edge.length)])
			
			# increment the poly.edges
			poly.edges <- poly.edges +1
			
			# move down to the new node
			current.polytomy <- new.node
			
			# get a new new.node
			new.node <- new.node - 1
		}
		
		# repair links below the new bifurcations
		phylo$edge[phylo$edge[,1]== orig.parent &
			phylo$edge[,2]== orig.polytomy] <- c(orig.parent, new.node +1)
	
	}
		
	root<-unique(setdiff(phylo$edge[,1],phylo$edge[,2]))
	phylo2<-phylo
	for (i in 1:length(phylo$edge[,1]))
	{
		for (j in 1:2)
		{
			phylo2$edge[i,j]<-if (phylo$edge[i,j]==root) "-1" else if (phylo$edge[i,j]=="-1") root else phylo$edge[i,j]
		}
	}
	
	write.tree(phylo2,file="temp.phy",multi.line=FALSE) #Modification to correct phylo structure
	phylo2<-read.tree(file="temp.phy") #Modification to correct phylo structure

	attr(phylo2, "force.binary") <- TRUE
	return(phylo2)
  }	
  }
}
\keyword{ manip}
