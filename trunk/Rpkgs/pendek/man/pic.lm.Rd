\name{pic.lm}
\alias{pic.lm}

\title{ Perform independent contrasts analysis with polytomous phylogenies }
\description{
Extends functionality of \code{pic} in five ways. First, it arranges for the phylogeny and 
data set to contain the same set of taxa. Second, it generates contrasts in all the 
variables requested. Third, it does not require a fully bifurcating phylogeny. Fourth, it 
optionally transforms branch lengths to optimise the evolutionary model assumed for each trait, 
reducing heteroscedasticity. Fifth, it optionally checks whether the fitted model is robust to 
the deletion of points having suspiciously high influence. }
\usage{
pic.lm(phylogeny, dataset, response, predictors, check.robust = FALSE, cutoff = 3, use.robust = FALSE, filters = NULL, kappa = NULL, file = NULL, plot=FALSE)
}

\arguments{
  \item{phylogeny}{An object of class \code{"phylo"}, which can contain polytomies}
  \item{dataset}{A data frame whose first column holds taxon names }
  \item{response}{Name of the column holding the response variable in the model to be fitted }
  \item{predictors}{ Vector of names of columns holding explanatory variables for the model}
  \item{check.robust}{ If \code{TRUE}, re-fits the model excluding points with studentised 
  residuals whose absolute value exceeds \code{cutoff} }
  \item{cutoff}{Number indicating how large a studentised residual must be before a point is excluded }
  \item{use.robust}{If \code{TRUE}, the returned model is that found with \code{check.robust}; by 
  default, the model returned is based on all points}
  \item{filters}{Vector of names holding the names of further columns holding dummy variables to permit 
  systematic exclusion of cases: taxa for which these columns contain missing values will be excluded from 
  the analysis. }
  \item{kappa}{Optional branch length transformation. If \code{kappa} is specified, branch lengths are 
  raised to the power of \code{kappa}. If \code{kappa} is not specified, it is optimised for each 
  variable in turn. }
  \item{file}{Optional file for contrasts to be written to. }
  \item{plot}{Is a model-checking plot required? Default = FALSE.}
}

\details{
  Polytomies are resolved arbitrarily into bifurcations separated by branches of zero length. The 
  contrasts computed at each of the \emph{n} bifurcations derived from a polytomy are given a weight 
  of 1/\emph{n} in the resulting analyses. This removes that arbitrary choice of a test variable 
  required by the CAIC program, while retaining validity in the face of incompletely resolved 
  phylogenies.
}

}
\value{
An object of class \code{lm}}

\references{Felsenstein, J. 1985  \emph{"American Naturalist"} \bold{"125"}, 1--15 for basic 
method of independent contrasts; Garland, T. Jr., Harvey, P. H., and Ives, A.R. 1992 \emph{Systematic 
Biology} \bold{41}, 18--32 for branch length optimization; Jones, K. E. and Purvis, A. 1997 
\emph{Functional Ecology} \bold{11},751--756 for treatment of influential observations}

\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\note{The algorithm for dealing with polytomies differs from that in the CAIC program, so results 
will likewise differ, even if \code{kappa} is fixed at 1}

\section{Warning }{At present there is limited error checking.}

\section{Warning }{Optimisation of \code{kappa} renders slopes more or less uninterpetable, but 
improves hypothesis-testing (reduces Type I and Type II error rates)}

\seealso{ \code{\link{test.linear}}, \code{\link{picmodsimp}}}

## The function is currently defined as
function(phylogeny, dataset, response, predictors, check.robust=FALSE, cutoff=3, use.robust=FALSE, filters=NULL, kappa=NULL,file=NULL, plot=FALSE){
#Written by Andy Purvis
#Models relationship between response variable and any  number of predictors.
#Optionally, other dummy variables can be used as filters to select subsets of the species
#If kappa is NULL, each variable has optimally transformed branch lengths.
#Otherwise, the specified kappa is used. Kappa = 0 gives equal branch lengths; kappa = 1
#is a null transformation (branch lengths remain unchanged)

#If check.robust = TRUE, repeat regression deleting contrasts whose studentised resisuals exceed ± cutoff

#response is a variable name; predictors and filters are vectors of variable names.

#requires ape
#requires MASS

#Change 3-3-06: Bug fix to make check.robust and use.robust work when phylogeny is fully binary
#Fix works by adding a weights column to the contrasts structure even when tree is binary

#Change 18-5-06: Bug fix to stop if there are not enough data to fit the model.
#Change 18-5-06: Bug fix to correct remaining problems caused by fully binary phylogeny.

#Change 19-5-06: The model-checking plot is now optional, and is not done by default.
#Change 19-5-06: Now checks the phylogeny is of class "phylo" and has an edge.length structure.

#Checks on phylogeny: added 19-5-06
if (class(phylogeny) != "phylo") 
    stop("object \"phy\" is not of class \"phylo\"")
if (is.null(phylogeny$edge.length)) 
    stop("your tree has no branch lengths: you may consider setting them equal to one, or using the function `compute.brlen'.")

variable.list<-c(response,predictors)
vars.and.filters<-c(response,predictors,filters)
dataset<-strip.spp.not.in.tree(dataset,phylogeny)
dataset<-strip.missing.data(dataset,vars.and.filters)
phylogeny<-strip.phylo(phylogeny,dataset)
tree.is.binary<-TRUE #Flag to keep track of whether input tree was binary

#check tree is big enough to go on with
nnodes<-length(phylogeny$edge[,1])-length(phylogeny$tip.label)+1
min.size<-length(predictors)+1

if (nnodes < min.size) stop("Not enough nodes for complexity of model")

if (!is.binary.tree(phylogeny)) 
{
	#Tree contains soft polytomies so must be made binary
	tree.is.binary<-FALSE
	poly.phylogeny<-phylogeny
	phylogeny<-force.binary.tree(poly.phylogeny)
	mapping<-map.nodes(phylogeny,poly.phylogeny)
}


nc<-length(dataset[,1])-1
contrasts<-array(data=NA,dim=c(nc,length(variable.list)))

if (is.null(kappa))
{
	#Need to optimise kappa for each variable in turn
	for (i in 1:length(variable.list))
	{
		contrasts[,i]<-optimal.pic(phylogeny, dataset,variable.list[i])
	}
} else {
	#Kappa is a specified constant, so transform branch lengths
	#But don't transform branches with no equivalent in the polytomous tree
	#The caveat is necessary if kappa = 0, when zero length branches become length 1
	# - and we only want that to happen for 'real zero-length branches
	
	phylogeny$edge.length<-phylogeny$edge.length^kappa
	#Now find 'artificial' branches coming from arbitrary resolution
	if(!tree.is.binary)
	{
		false.nodes<-subset(mapping[,2],mapping[,4]==FALSE)
		indices<-which(is.element(phylogeny$edge[,2],false.nodes))
		#Set these branches to length zero
		phylogeny$edge.length[indices]<-0
	}
	
	for (i in 1:length(variable.list))
	{
		species.data<-dataset[,names(dataset)==variable.list[i]]
		names(species.data)<-dataset[,1]
		contrasts[,i]<-pic(species.data,phylogeny)
	}
}

contrasts<-as.data.frame(contrasts,row.names=attributes(contrasts[,1],which="names"))
names(contrasts)<-variable.list

contrasts$node.in.original.tree<-rep(NA,length(contrasts[,1])) #Added 3-3-06 to harmonise structure
contrasts$regression.weights<-rep(1,length(contrasts[,1])) #Added 3-3-06 to harmonise structure

if (tree.is.binary==FALSE) 
{
	#Tree contains soft polytomies so add weights and mapping
	contrasts$node.in.original.tree<-mapping[,1]
	contrasts$regression.weights<-as.numeric(mapping[,3])
}


if(!is.null(file))					#User has requested output file and specified name
	{
	write.table(file=file, contrasts, quote=F, sep="\t", append=F, row.names=T, col.names=T) #Save results
	}


predictors<-paste("contrasts$",predictors,sep="")
response<-paste("contrasts$",response,sep="")
tmp<-paste(predictors,collapse=" + ")
fmla<-paste(response,"~",tmp,"- 1")
print(fmla)

if (tree.is.binary==FALSE)
{
	model1<-lm(as.formula(fmla),weights=contrasts$regression.weights)
	model1$df.residual<-round(sum(contrasts$regression.weights)-length(predictors))
	print(summary.pic.lm(model1))
} else {
	model1<-lm(as.formula(fmla),weights=contrasts$regression.weights)
	print(summary(model1))
}

if (plot==TRUE & check.robust==FALSE){
	#19/5/06: Now an option, used to be automatic
	par(mfrow=c(2,2))
	plot(model1)
	par(mfrow=c(1,1))
}



print("Any warnings?")
print(warnings())
if (check.robust==TRUE)
{
	# browser()
	if (length(model1$resid) - model1$rank <= 1) stop("Too few observations to compute studentized residuals") #Added 19-5-06
	
	print(paste("Checking effect of deleting contrasts with absolute studentised residuals >", cutoff,"..."))
	stud.res<-my.stud.res(model1)
	to.delete<-which(abs(stud.res)>cutoff)
	print(paste("Deleting",length(to.delete),"contrasts..."))
	to.keep<-rep(1,length(contrasts[,1]))
	to.keep[to.delete]<-0
	contrasts<-subset(contrasts,to.keep==1)
	
	if (sum(contrasts$regression.weights) < min.size) stop("Not enough observations remaining to fit model") #Added 18-5-06
	if (tree.is.binary==FALSE)
	{
		model1r<-lm(as.formula(fmla),weights=contrasts$regression.weights)
		model1r$df.residual<-round(sum(contrasts$regression.weights)-length(predictors))
		print(summary.pic.lm(model1r))
	} else {
		model1r<-lm(as.formula(fmla),weights=contrasts$regression.weights)
		print(summary(model1r))
	}
	
	if (plot==TRUE){
		#19-5-06: Now an option, used to be automatic
		par(mfrow=c(2,2)) 
		plot(model1r)
		par(mfrow=c(1,1))
		}
		
	if (use.robust==TRUE)
	{
		model1<-model1r
		if(!is.null(file))					#User has requested output file and specified name
		{
			write.table(file=file, contrasts, quote=F, sep="\t", append=F, row.names=T, col.names=T) #Overwrite with pruned contrasts file
		}

	}
}

return(model1)
}

\keyword{regression}
\keyword{htest}
\keyword{dplot}
\keyword{hplot}
\keyword(robust)