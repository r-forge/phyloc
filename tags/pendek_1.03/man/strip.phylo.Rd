\name{strip.phylo}
\alias{strip.phylo}

\title{ Remove terminals from a phylogeny that are not present in a data frame of comparative data }
\description{
To use \code{pic}, \code{pic.lm} and variants, the phylogeny and data set need to have matching sets 
of taxa. This function is one of several that achieves the match.}
\usage{
strip.phylo(phylogeny, dataset)
}

\arguments{
  \item{phylogeny}{An object of class \code{"phylo"} }
  \item{dataset}{A data frame whose first column holds taxon names }
}
\details{
An easy problem that might be encountered is that taxon names in the data set often have 
genus and species names separated by a space, whereas an underscore is the more likely 
separator in a phylogeny; using (e.g.) \code{gsub} on the data frame taxon names can fix 
this simply.}
\value{
An object of class \code{"phylo"}, with only those terminal taxa that are also 
present in column 1 of \code{dataset}}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking on the phylogeny passed to the function}

\seealso{\code{\link{strip.missing.data}}, \code{\link{strip.spp.not.in.tree}},\code{\link[ape]{pic}} }

## The function is currently defined as
function(phylogeny, dataset){

#Written by Andy Purvis
#Drops tips from a phylogeny that are not present in the data set

#requires ape
old<-sum(phylogeny$edge>0)
not.in.data<-setdiff(phylogeny$tip.label,dataset[,1])
if (length(not.in.data)>0) phylogeny<-drop.tip(phylogeny,not.in.data)
new<-sum(phylogeny$edge>0)
print(paste("Phylogeny stripped from",old,"to",new,"species."))

phylogeny
  }
}
\keyword{manip }
