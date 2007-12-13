"my.stud.res" <-
function(model){
#Written because of strange behaviour of MASS function stdres (and hence studres).
#See Venables and Ripley MASS book, pp 204-5

	std.res<-my.std.res(model)
	p<-length(coef(model))
	n<-length(model$residuals)
	stud.res<-std.res/sqrt((n-p-std.res^2)/(n-p-1))
	
	stud.res
}

