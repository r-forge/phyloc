"plot.hcd" <-
function(x, ...){

# the log = "y" call should be particularly useful so
# make sure it is passed to plot correctly.

if(is.null(x$max)) x$min <- x$max <- x$means

# get limits
if(is.null(x$real)){
	max.clade.size <- max(x$max)
	min.clade.size <- min(x$min)
	main <- paste("HCD of", x$species, "species among", x$taxa,"from", x$reps, "simulation replicates")	
	}
else{
	max.clade.size <- max(x$max,x$real)
	min.clade.size <- min(x$min,x$real)
	main <- paste("HCD for : ", x$data, ".\n",x$species," species among ", x$taxa," taxa from ", x$reps, " simulation replicates.",sep="")	
}


plot(x$means,ylim=c(min.clade.size,max.clade.size),type="l",main=main,xlab="Taxa",ylab="Number of species",...)

if(!is.null(x$real))
	lines(x$real,col="red")
if(!is.null(min)){			# assume max is also present
	lines(x$min,lty=3)
	lines(x$max,lty=3)
}

# put a legend on it?

}

