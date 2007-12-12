"clifford.acf" <-
function(x, ew.wrap=FALSE){

    # This is slow by comparison to the Fortran version:
    # ~ 45 minutes for a 360x152 matrix as compared to ~ 7 minutes

    # get the dimensions and the output matrices
    nr <- dim(x)[1]
    nc <- dim(x)[2]
    nok <- acf <- matrix(0, ncol=(nc*2)-1, nrow=(nr*2)-1)
        
    # create a 0/1 matrix showing non-NA values
    not.na.x <- !(is.na(x))
    # and a copy with NA set to zero to avoid n+NA -> 0
    xx <- x
    xx[is.na(xx)] <- 0
    
	# sum the products of each cell against each lag...
	for(r in 1:nr){
		for(c in 1:nc){
			#  ... as long as the cell isn't NA
			if(!is.na(x[r,c])){
				indr <- ((nr + 1) - r):(((nr + 1) - r) + (nr - 1))
				indc <- ((nc + 1) - c):(((nc + 1) - c) + (nc - 1))
				acf[indr,indc] <- acf[indr, indc] + (xx*xx[r,c])
				nok[indr,indc] <- nok[indr, indc] + not.na.x           
			}
		}
	}   
    
    if(ew.wrap){
    	
    	acf.nc <- ((nc*2)-1) 

		# find the number of columns to wrap from each side
        # arbitrarily breaking even numbers of columns
        cols.to.wrap <- acf.nc - nc
        left.wrap <- ceiling(cols.to.wrap/2) 
        right.wrap <- floor(cols.to.wrap/2)
        
        # calculate vectors of the old and new positions
        lo <- 1:left.wrap
        ro <- (acf.nc - right.wrap + 1):acf.nc
        ln <- (acf.nc - right.wrap + 1 - left.wrap):(acf.nc - right.wrap)
        rn <- (left.wrap + 1):(left.wrap + right.wrap)
        
        acf[,rn] <- acf[,rn] + acf[,ro]
		acf[,ln] <- acf[,ln] + acf[,lo]
		acf[,lo] <- 0
		acf[,ro] <- 0
		
        nok[,rn] <- nok[,rn] + nok[,ro]
		nok[,ln] <- nok[,ln] + nok[,lo]
		nok[,lo] <- 0
		nok[,ro] <- 0

       
    }
    
    acf <- acf/nok
    
    # need to allow for lags that have no representations in
    # the dataset - the averaging above will do 0/0 and hence NaN
    acf[is.nan(acf)]  <- 0
    
    return(list(acf=acf, nok=nok))
}

