"print.hcd" <-
function(x, ...){

cat("\nHollow Curve Distribution of species among higher taxa\n")
if(!is.null(x$real))
	cat("        HCD parameters from data:",x$data)
cat("\n\nNumber of species:",x$species,"\n")
cat("Number of higher taxa:", x$taxa,"\n")
cat("Number of replicates:", x$reps,"\n\n")

if(!is.null(x$hi)){
	if(x$num.hi == 0){
		hiprob <- paste("<",1/x$reps)
	} else {
		hiprob <- formatC(x$num.hi/x$reps, digits=2)
	}
	cat("Probability of largest sub-taxon (", x$hi, "): ", hiprob, "\n", sep="")
}

if(!is.null(x$lo)){
	if(x$num.lo == 0){
		loprob <- paste("<",1/x$reps)
	} else {
		loprob <- formatC(x$num.lo/x$reps, digits=2)
	}
	cat("Probability of smallest sub-taxon (", x$lo, "): ", loprob, "\n", sep="")
}

cat("\n")

}

