\name{strip.spp.not.in.tree}
\alias{strip.spp.not.in.tree}

\title{Remove rows from a data frame of comparative data those rows corresponding to taxa 
not in the phylogeny }
\description{
To use \code{pic}, \code{pic.lm} and variants, the phylogeny and data set need to have matching sets 
of taxa. This function is one of several that achieves the match.}
}
\usage{
strip.spp.not.in.tree(dataset, phylogeny)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{"phylo"} }
  \item{dataset}{A data frame whose first column holds taxon names }
}
\details{
An easy problem that might be encountered is that taxon names in the data set often have 
genus and species names separated by a space, whereas an underscore is the more likely 
separator in a phylogeny; using (e.g.) \code{gsub} on the data frame taxon names can fix 
this simply.}
}
\value{
A \code{data.frame}, with only those rows whose taxa are also in \code{phylogeny}}


\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking on the phylogeny passed to the function}

\seealso{\code{\link{strip.missing.data}}, \code{\link{strip.phylo}},\code{\link[ape]{pic}} }

## The function is currently defined as
function(dataset, phylogeny){

#Written by Andy Purvis
#Omits species from a data frame (dataset) that are not present in a phylogeny
#Presumes (as do other functions here) that species names are in the first column of dataset

#requires ape
old<-length(dataset[,1])
matches<-match(dataset[,1],phylogeny$tip.label,nomatch=0)
dataset<-subset(dataset,matches>0)
new<-length(dataset[,1])
print(paste("Dataset stripped from",old,"to",new,"species."))

dataset
  }
}
\keyword{ manip}