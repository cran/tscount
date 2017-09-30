vcov.tsglm <- function(object, ...){
  if(is.null(object$info.matrix)) stop("No information matrix provided. Argument 'object' must be the output of a call\nto the function 'tsglm' with argument 'info' not set to \"none\"")
  invertedinfo <- invertinfo(object$info.matrix, stopOnError=TRUE)$vcov
  if(is.null(object$info.matrix_corrected)){
    warning("No corrected information matrix provided, so that the information matrix of a\nPoisson model is used. If a negative binomial model was fitted, argument 'object'\nmust be the output of a call to the function 'tsglm' with argument 'info' set\nto \"score\"")
    result <- invertedinfo
  }else{
    result <- invertedinfo %*% object$info.matrix_corrected %*% invertedinfo #sandwich-type formula (is equal to invertedinfo in case of a Poisson distribution)
    #result <- (result + t(result))/2 #ensures that the resulting matrix is exactly symmetric, because the previous matrix multiplication may result in numerical deviations from its symmetry
  }
  return(result)
}
