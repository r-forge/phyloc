"export.hcd" <-
function(x, outfile){

if(class(x) != "hcd") stop("Requires an object of class 'hcd'")

outframe <- data.frame(means = x$means)

if(!is.null(x$real))
	outframe <- cbind(outframe, real=x$real)
	
if(!is.null(x$min))
	outframe <- cbind(outframe, min=x$min)
	
if(!is.null(x$max))
	outframe <- cbind(outframe, max=x$max)

con <- file(outfile, open="w")

cat("# Hollow Curve Distribution using the broken stick model\n# ", 
	x$reps , " simulations of ", x$species ," species among ", 
	x$taxa , " taxa\n", sep="",file=con)
	
if(!is.null(x$real))
	cat("# Based on data from ", x$data, "\n", sep="",file=con, append=TRUE)
	
if(!is.null(x$hi)){
	if(x$num.hi == 0){
		hiprob <- paste("<",1/x$reps)
	} else {
		hiprob <- formatC(x$num.hi/x$reps, digits=2)
	}
	cat("# Probability of largest sub-taxon (", x$hi, "): ", hiprob, "\n", sep="", file=con, append=TRUE)
}

if(!is.null(x$lo)){
	if(x$num.lo == 0){
		loprob <- paste("<",1/x$reps)
	} else {
		loprob <- formatC(x$num.lo/x$reps, digits=2)
	}
	cat("# Probability of smallest sub-taxon (", x$lo, "): ", loprob, "\n", sep="", file=con, append=TRUE)
}

cat(paste(names(outframe), collapse="\t"), "\n", sep="", append=TRUE, file=con)
close(con)

write.table(outframe, file=outfile, row.names=FALSE, append=TRUE, col.names=FALSE)

}

