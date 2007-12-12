
#This is the R function to call the ReadWithNCL C++ function
#other params can be added to the argument list, and then
#added to the params list below
ReadWithNCL <- function(fileToRead){

# Most of the input parameter checking here is not really
# necessary because it is done in the Rcpp code.
  
  params <- list(filename=fileToRead)
  
  # Check that params is properly formatted.
  if(!is.list(params) || length(params) == 0) {
    stop("The params parameter must be a non-empty list");
  }

  # Finally ready to make the call...
    val <- .Call("ReadWithNCL", params,
               PACKAGE="ReadWithNCL")

  # Define a class for the return value so we can control what gets
  # printed when a variable assigned this value is typed on a line by itself.
  # This has the effect of calling the function print.RcppExample(). The
  # function (defined below) simply prints the names of the fields that are
  # available. Access each field with val$name.
  class(val) <- "ReadWithNCL"
  
  val
}


RcppExample <- function(params, nlist, numvec, nummat, df, datevec, stringvec,
                        fnvec, fnlist) {

# Most of the input parameter checking here is not really
# necessary because it is done in the Rcpp code.
  
  # Check that params is properly formatted.
  if(!is.list(params) || length(params) == 0) {
    stop("The params parameter must be a non-empty list");
  }

  # Check nlist
  if(!is.list(nlist) || length(nlist) == 0) {
    stop("The nlist parameter must be a non-empty list");
  }
  if(length(nlist) != length(names(nlist))) {
    stop("The values in nlist must be named")
  }
  if(!is.numeric(unlist(nlist))) {
    stop("The values in nlist must be numeric")
  }

  # Check numvec argument
  if(!is.vector(numvec)) {
    stop("numvec must be a vector");
  }

  # Check nummat argument
  if(!is.matrix(nummat)) {
    stop("nummat must be a matrix");
  }
  
  # Finally ready to make the call...
  val <- .Call("Rcpp_Example", params, nlist, numvec, nummat,
               df, datevec, stringvec, fnvec, fnlist,
               PACKAGE="RcppTemplate")

  # Define a class for the return value so we can control what gets
  # printed when a variable assigned this value is typed on a line by itself.
  # This has the effect of calling the function print.RcppExample(). The
  # function (defined below) simply prints the names of the fields that are
  # available. Access each field with val$name.
  class(val) <- "RcppExample"
  
  val
}

print.RcppExample <- function(x,...) {
  cat('Names defined in RcppExample return list:\n')
  cat('(Use result$name to access)\n')
  namevec <- names(x)
  for(i in 1:length(namevec))
    cat(i, ': ', namevec[i], '\n')
}

