"clademat2vcv" <-
function(clmat, method="loop"){

	n.tips <- dim(clmat$clade.matrix)[2]
	V <- matrix(NA, n.tips, n.tips)
	
	
	switch(method,
	"loop" = for(i in 1:n.tips){
					for(j in 1:i){		
						V[i,j] <- sum(clmat$edge.length[clmat$clade.matrix[,i] & clmat$clade.matrix[,j]])
					}
				},
	"bycolumn" =  for(i in 1:n.tips){
					V[i,] <- colSums(clmat$edge.length * (clmat$clade.matrix[,i] & clmat$clade.matrix))
				  })
	
	dimnames(V) <- list(clmat$tip.label, clmat$tip.label)
	return(V)

}

