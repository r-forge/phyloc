makeNamedVec <- function(x, colName) {
	vec <- x[,colName]
	names(vec) <- row.names(x)
	vec
}
