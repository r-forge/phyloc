"my.std.res" <-
function(model){
#Generates standardised residuals for a model
#Written because of strange behaviour of MASS function stdres.
#See Venables and Ripley MASS book, p 204

	hats<-lm.influence(model)$hat
	s<-summary.pic.lm(model)$sigma
	std.res<-(model$residuals*sqrt(model$weights))/(s*sqrt(1-hats))
	
	std.res
}

