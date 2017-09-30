simcoefs <- function(...) UseMethod("simcoefs")

simcoefs.tsglm <- function(fit, method=c("bootstrap", "normapprox"), B=1, parallel=FALSE, ...){
  stopifnot(
    length(B)==1,
    B%%1==0,
    B>=1
  )
  method <- match.arg(method)
  if(method=="bootstrap"){
    simfit <- function(seed, fit, ...){
      set.seed(seed)
      ts_sim <- tsglm.sim(fit=fit)$ts
      suppressWarnings(fit_sim <- tsglm(ts=ts_sim, model=fit$model, xreg=fit$xreg, link=fit$link, distr=fit$distr, score=FALSE, info="none", ...))
      if(fit$distr=="nbinom" && fit_sim$distr=="poisson") fit_sim$distrcoefs <- c(size=NA)
      result <- c(coef(fit_sim), sigmasq=fit_sim$sigmasq, fit_sim$distrcoefs)
      return(result)
    }
    seeds <- sample(1e+9, size=B)
    if(parallel){
      Sapply <- function(X, FUN, ...) parSapply(cl=NULL, X=X, FUN=FUN, ...)
    }else{
      Sapply <- sapply
    }
    coefs <- t(Sapply(seeds, simfit, fit=fit, ...))
    result <- list(coefs=coefs)
  }
  if(method=="normapprox"){
    rmvnorm_stable <- function(n, mean=rep(0, nrow(sigma)), sigma=diag(length(mean))){
      #Function for stable generation of random values from a multivariate normal distribution (is robust against numerical deviations from symmetry of the covariance matrix. Code is taken from function rmvnorm in the package mvtnorm and modified accordingly.
      if(length(mean) != nrow(sigma)) stop("mean and sigma have non-conforming size")
      ev <- eigen(sigma, symmetric=TRUE)
      ev$values[ev$values < 0] <- 0
      R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
      centred <- matrix(rnorm(n=n*ncol(sigma)), nrow=n, byrow=TRUE) %*% R
      result <- sweep(centred, 2, mean, "+")
      colnames(result) <- names(mean)
      return(result)
    }
    f <- 1.1 #one could choose this factor according to the probability of a parameter from the multivariate normal distribution to be outside the parameter space
    coefs <- rmvnorm_stable(n=ceiling(f*B), mean=coef(fit), sigma=vcov(fit))
    repeat{
      valid_coefs <- apply(coefs, 1, function(x) tsglm.parametercheck(tsglm.parameterlist(paramvec=x, model=fit$model), link=fit$link, stopOnError=FALSE))
      if(sum(valid_coefs) >= B) break
      coefs <- rbind(coefs, rmvnorm_stable(n=ceiling((B-sum(valid_coefs))*f/mean(valid_coefs)), mean=coef(fit), sigma=vcov(fit)))
    }
    use_coefs <- which(valid_coefs)[1:B]
    coefs <- coefs[use_coefs, , drop=FALSE]
    n_invalid <- max(use_coefs) - B
    distrcoefs_matrix <- matrix(rep(c(sigmasq=fit$sigmasq, fit$distrcoefs), B), byrow=TRUE, nrow=B)
    colnames(distrcoefs_matrix) <- c("sigmasq", names(fit$distrcoefs))
    coefs <- cbind(coefs, distrcoefs_matrix)  
    result <- list(coefs=coefs, n_invalid=n_invalid)
  }
  return(result)
}        
