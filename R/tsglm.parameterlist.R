tsglm.parameterlist <- function(paramvec, model){
  p <- length(model$past_obs)
  P <- seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  q <- length(model$past_mean)
  Q <- seq(along=numeric(q)) #sequence 1:p if p>0 and NULL otherwise
  r <- length(paramvec) - (1+p+q)
  R <- seq(along=numeric(r))
  names(paramvec) <- NULL 
  result <- list(intercept=paramvec[1], past_obs=paramvec[1+P], past_mean=paramvec[1+p+Q], xreg=paramvec[1+p+q+R])
  return(result)
}
  