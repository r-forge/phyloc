makeNamedVec <- function(x, colnum) {
	vec <- x[,colnum]
	names(vec) <- row.names(x)
}
