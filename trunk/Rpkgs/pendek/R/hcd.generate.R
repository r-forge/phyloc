"hcd.generate" <-
function(species, taxa, reps = 1000, hi = NULL , lo = NULL, minmax = FALSE){

# approx. 500 cycles a second on Mac G3 266MHz 160 Mb RAM

	reps.store <- array(data=0, dim=c(taxa, reps))
	
	if(minmax){
		min.store <- max.store <- vector("numeric",length = taxa)
		min.store <- species #intialized so vector contains large numbers
	}

	num.hi <- 0
	num.lo <- 0

	for(replicates in 1:reps){

		breakpoints <- sort(sample(1:(species-1), taxa-1))
		breaklengths <- sort(diff(c(0,breakpoints,species)), decreasing=TRUE)
		# in decreasing size order
				
		reps.store[,replicates] <- breaklengths
		
		if(minmax){
			max.store <- pmax(max.store,breaklengths)
			min.store <- pmin(min.store,breaklengths)
		}
		
		# calculate probs as proportion of simulated taxa that 
		# exceed criteria rather than proportion of simulated runs?
		if(!is.null(hi))
			if(any(breaklengths>=hi)) num.hi <- num.hi + 1
		if(!is.null(lo))
			if(any(breaklengths<=lo)) num.lo <- num.lo + 1

	}

means <- apply(reps.store, 1, sum)/reps

OUTPUT <- list(species = species, taxa = taxa, reps = reps, means = means)

if(minmax){
	OUTPUT$min <-  min.store
	OUTPUT$max <-  max.store
	}
if(!is.null(hi))
	OUTPUT <- c(OUTPUT, hi = hi, num.hi = num.hi)
if(!is.null(lo))
	OUTPUT <- c(OUTPUT, lo = lo, num.lo = num.lo)
class(OUTPUT)<-"hcd"
return(OUTPUT)

}

