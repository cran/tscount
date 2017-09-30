#Density function of the conditional distribution:
ddistr <- function(x, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, ...){
  distr <- match.arg(distr)
  result <- switch(distr,
    "poisson"=dpois(x, lambda=meanvalue, ...),  
    "nbinom"=dnbinom(x, mu=meanvalue, size=distrcoefs[[1]], ...)
  )  
  return(result)
}

#Cumulative distribution function of the conditional distribution:
pdistr <- function(q, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, ...){
  distr <- match.arg(distr)
  result <- switch(distr,
    "poisson"=ppois(q, lambda=meanvalue, ...),  
    "nbinom"=pnbinom(q, mu=meanvalue, size=distrcoefs[[1]], ...)
  )  
  return(result)
}

#Quantile function of the conditional distribution:
qdistr <- function(p, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, ...){
  distr <- match.arg(distr)
  result <- switch(distr,
    "poisson"=qpois(p, lambda=meanvalue, ...),  
    "nbinom"=qnbinom(p, mu=meanvalue, size=distrcoefs[[1]], ...)
  )  
  return(result)
}

#Random number generation from conditional distribution:
rdistr <- function(n, meanvalue, distr=c("poisson", "nbinom"), distrcoefs){
  distr <- match.arg(distr)
  result <- switch(distr,
    "poisson"=rpois(n, lambda=meanvalue),  
    "nbinom"=rnbinom(n, mu=meanvalue, size=distrcoefs[[1]])
  )  
  return(result)
}

#Standard deviation of the conditional distribution:
sddistr <- function(meanvalue, distr=c("poisson", "nbinom"), distrcoefs){
  distr <- match.arg(distr)
  result <- switch(distr,
    "poisson"=sqrt(meanvalue),  
    "nbinom"=sqrt(meanvalue + meanvalue^2/distrcoefs[[1]])
  )  
  return(result)
}

#Anscombe residuals:
ardistr <- function(response, meanvalue, distr=c("poisson", "nbinom"), distrcoefs){
  result <- switch(distr,
    "poisson"=3/2*(response^(2/3)-meanvalue^(2/3))/meanvalue^(1/6),  
    "nbinom"=(3/distrcoefs[["size"]]*((1+response*distrcoefs[[1]])^(2/3) - (1+meanvalue*distrcoefs[[1]])^(2/3)) + 3*(response^(2/3)-meanvalue^(2/3))) / (2*(meanvalue+meanvalue^2*distrcoefs[[1]])^(1/6))
  )  
  return(result)
}

checkdistr <- function(distr=c("poisson", "nbinom"), distrcoefs){
  distr <- match.arg(distr)
  if(distr=="nbinom"){
    if(missing(distrcoefs) || length(distrcoefs)!=1) stop("For the negative binomial parameter (only) the dispersion parameter 'size' has to be provided in argument 'distrcoefs'")
    if(distrcoefs[[1]]<=0) stop("The additional dispersion parameter for the negative binomial distribution has to be greater than zero")
  }
}
