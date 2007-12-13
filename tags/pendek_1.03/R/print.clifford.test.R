"print.clifford.test" <-
function(x, ...){

cat("\nCorrelation test accounting for autocorrelation (Clifford et al., 1989)\n\nData - Matrix A:", x$A$name, 
    "\n     - Matrix B:", x$B$name, "\n\n")
print(as.data.frame(x[4:12], row.names=""))
cat("\n")
}

