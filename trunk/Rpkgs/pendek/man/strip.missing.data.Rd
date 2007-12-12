\name{strip.missing.data}
\alias{strip.missing.data}

\title{ Remove incomplete cases from a comparative data set}
\description{
To use \code{pic}, \code{pic.lm} and variants, the phylogeny and data set need to have matching sets 
of taxa. This function is one of several that achieves the match.}
}
\usage{
strip.missing.data(dataset, variable.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{A data frame whose first column holds taxon names  }
  \item{variable.list}{ A vector holding the \emph{names} of columns to be used in a particular analysis }
}
\details{
Taxa can only be used in independent contrasts analyses if their values of /emph{all} the traits 
specified in the model are in the data set.  This function prunes out any incomplete cases.}
\value{
A \code{data.frame}, with only those rows for which no relevant variable has a missing value}
}


\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking on the arguments passed to the function}

\seealso{\code{\link{strip.spp.not.in.tree}}, \code{\link{strip.phylo}},\code{\link[ape]{pic}} }

## The function is currently defined as
function(dataset, variable.list){

#Written by Andy Purvis
#Omits species from a dataset that have missing values for any of the variables in variable.list
#variable.list is a vector of variable names of interest

old<-length(dataset[,1])
for (i in 1:length(variable.list))
{
	dataset<-subset(dataset,!is.na(dataset[names(dataset)==variable.list[i]]))
  }
new<-length(dataset[,1])
print(paste("Dataset stripped from",old,"to",new,"species."))

dataset
  }
}
\keyword{manip}
