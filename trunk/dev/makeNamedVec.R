makeNamedVec <- function(x, colID) {
	vec <- x[,colID]
	names(vec) <- row.names(x)
	vec
}
