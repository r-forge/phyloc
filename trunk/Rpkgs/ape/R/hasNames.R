hasNames <- function(x) {

if (is.vector(x)) {
	if (length(names(x)) == 0) return(FALSE)
	else return(TRUE)
}
else if (is.data.frame(x) || is.matrix(x)) {
	if (length(row.names(x)) == 0) return(FALSE)
	else return(TRUE)
}

}
