tsglm.distrfit <- function(object, distr){ 
  if(distr=="poisson"){
    distrcoefs <- NULL
    sigmasq <- 0
  }
  if(distr=="nbinom"){  
    fitval <- object$fitted.values
    ts <- object$response
    n <- object$n_eff
    m <- length(object$coefficients)  
    #Pearson type estimator:
    find_root <- function(v) sum((ts-fitval)^2/(fitval*(1+fitval/v))) - n + m
    root <- try(uniroot(f=find_root, interval=c(0, 1e100)), silent=TRUE)  
    if(class(root)=="try-error"){
      distr <- "poisson"
      distrcoefs <- NULL
      sigmasq <- 0
      warning("The dispersion parameter of the negative binomial distribution cannot be\nestimated. This indicates that there is no or only very little overdispersion\nin the data. The Poisson distribution with argument 'distr' set to \"poisson\"\nwas fitted instead.")
    }else{
      distrcoefs <- c(size=root$root)
      sigmasq <- 1/root$root
    }
  }
  result <- list(distr=distr, distrcoefs=distrcoefs, sigmasq=sigmasq)
  return(result)
}
