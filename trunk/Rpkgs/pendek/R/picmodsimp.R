"picmodsimp" <-
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

