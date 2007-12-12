\name{picmodsimp}
\alias{picmodsimp}

\title{Automated model simplification for independent contrasts multiple regressions}
\description{
First a multiple regression model using contrasts in all \code{predictors}. Then drops the least 
significant predictor variables, with higher-order terms being dropped before main effects, until 
all terms remaining in the model are either significant or are present in significant higher-order 
terms. The algorithm used is greedy in the presence of missing values; see details.}
\usage{
picmodsimp(phylogeny, dataset, response, predictors, main.effects, p.crit = 0.05, check.robust = FALSE, cutoff = 3, use.robust = FALSE, filters = NULL, kappa = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{An object of class \code{"phylo"}, which can contain polytomies}
  \item{dataset}{A data frame whose first column holds taxon names, second column holds the variable specified in \code{response}, 
  with the main effects occupying the next columns. Names of higher-order variables must exactly incorporate the name of their 
  corresponding main effect(s).}
  \item{response}{Name of the column holding the response variable in the model to be fitted; NOTE - must be the first column 
  after the taxon names. }
  \item{predictors}{ Vector of names of columns holding explanatory variables for the model}
  \item{main.effects}{ numeric: number of main effects in starting set of variables. These variables must occupy the columns 
  immediately after the response variable in the data set. }
  \item{p.crit}{ numeric: critical p-value for keeping terms in the model during simplification. }
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
}
\details{
The data set must be set up precisely as specified. The procedure is greedy in the face of missing values in the data, and 
does not guarantee finding the best minimum adequate model. The greed can be ameliorated by using \code{with.fixed.predictors} 
to re-enter all dropped variables in turn into the model returned by this function, to see if any are significant.}
\value{
An object of class \code{lm}}
}
\references{ Purvis, A., Gittleman, J. L., Cowlishaw, G. & Mace, G. M. 2000 \emph{Proc. R. Soc. Lond. B} \bold{267}, 1947--1952 for search approach }
\author{Andy Purvis}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\section{Warning }{At present there is no error checking}

\seealso{\code{\link{pic.lm}}, \code{\link{with.fixed.predictors}} }

## The function is currently defined as
function(phylogeny, dataset, response, predictors, main.effects, p.crit=0.05,check.robust=FALSE,cutoff=3,use.robust=FALSE,filters=NULL,kappa=NULL){
#Fits linear model on contrasts in all variables in response ~ predictors.
#Then drops least useful term, recomputes contrasts if necessary, and so on until a MAM is found.

#User specifies data frame, the name of the response variable (which is a bit unnecessary, since it is constrained
#to be the first column, and how many main effects there are (these must occupy the columns immediately after the
#response variable). Subsequent columns hold interactions, higher powers etc; THE NAMES MUST INCORPORATE EXACTLY THE
#NAMES OF COMPONENT MAIN EFFECTS!

#Each main effect is attributed the lowest p-value associated with any term containing it, so that non-main-effects 
#are dropped in preference to main effects.

#The term with the highest p-value attributed to it is dropped (in event of a tie, the term later in the list is
#dropped, again to keep main effects as long as possible) until all terms have significant p-values attributed to them.

#The user can override the default critical p-value of 0.05 by specifying another p.crit in the call

#requires ape

model<-pic.lm(phylogeny,dataset,response,predictors,check.robust=check.robust,cutoff=cutoff,use.robust=use.robust,kappa=kappa,filters=filters)

still.in<-predictors
main<-c(rep(TRUE,main.effects),rep(FALSE,length(still.in)-main.effects))
p.value<-summary.pic.lm(model)$coefficients[,4]

max.p<-max(p.value)
while (max.p > p.crit)
{

	for (i in 1:main.effects)
	{
		relevant.terms<-which(regexpr(still.in[i],still.in)!=-1)
		relevant.p<-p.value[relevant.terms]
		min.p<-min(relevant.p)
		p.value[i]<-min.p #Attribute to main effect the lowest p-value of any term including it
	}
	
	to.lose<-max(which(p.value==max(p.value))) #In event of tie, take variable that is later in variable list
	print(paste("Least useful term remaining:",still.in[to.lose],"; p.value =",max(p.value)))
	print("")
	
	if (max(p.value)>p.crit)
	{
		to.lose.name<-still.in[to.lose]
		print(to.lose.name)
		still.in<-still.in[-to.lose]

		if (to.lose>main.effects)
		{
			#Don't need to recompute contrasts, as sample size can't have changed
			model <- update(model, as.formula(paste(".~. -", to.lose.name)))
		}
		if (to.lose<=main.effects)
		{
			#Need to recompute contrasts as sample size could have changed
			main.effects<-main.effects-1
			print("**************************")
			print(paste("**************************   LOSING A MAIN EFFECT - NOW HAVE ",main.effects))
			print("**************************")
	
			model<-pic.lm(phylogeny,dataset,response,still.in,check.robust=check.robust,cutoff=cutoff,use.robust=use.robust,kappa=kappa,filters=filters)
		}
		p.value<-summary.pic.lm(model)$coefficients[,4]
		max.p<-max(p.value)
	}
	if (max(p.value)<=p.crit) 
	{
		print("****MODEL IS MINIMUM-ADEQUATE")
		max.p<-max(p.value)
	}
  }
print(summary.pic.lm(model))

model
  }
}
\keyword{ regression}
\keyword{ robust}