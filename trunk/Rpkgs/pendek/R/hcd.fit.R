"hcd.fit" <-
function(x,...){

if(is.vector(x)&&is.numeric(x)){
	species <- sum(x)
	taxa <- length(x)
	hi <- max(x)
	lo <- min(x)
	fitted.hcd <- hcd.generate(species,taxa, hi = hi, lo = lo,...)
	fitted.hcd$real <- sort(x, decreasing=TRUE)
	fitted.hcd$data <- deparse(substitute(x))
	}
else{
	stop("Requires a vector of integers.")
	}
	
return(fitted.hcd)
}

