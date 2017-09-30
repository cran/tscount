se <- function(object, ...) UseMethod("se")

se.tsglm <- function(object, B, parallel=FALSE, level=0.95, ...){
  tsglm.check(object)
  stopifnot(
    length(level)==1,
    !is.na(level),
    level<1,
    level>0
  )
  a <- (1-level)/2
  est <- c(coef(object), sigmasq=if(object$distr=="poisson") NULL else object$sigmasq)
  if(missing(B)){
    covmatrix <- vcov(object)
    variances <- diag(covmatrix)
    stderrors <- c(sqrt(variances), sigmasq=if(object$distr=="poisson") NULL else NA)
    ci <- cbind(lower=est-qnorm(1-a)*stderrors, upper=est+qnorm(1-a)*stderrors)
    result <- list(est=est, se=stderrors, ci=ci, level=level, type="normapprox")
  }else{
    bootstrap_coefs <- simcoefs(object, method="bootstrap", B=B, parallel=parallel, ...)$coefs[, seq(along=est), drop=FALSE]
    if(object$distr!="poisson"){
      n_invalid <- sum(bootstrap_coefs[, "sigmasq"]==0)
      if(n_invalid>0) warning(paste("The overdispersion coefficient 'sigmasq' could not be estimated\nin", n_invalid, "of the", B, "replications. It is set to zero for these\nreplications. This might to some extent result in a biased estimation\nof its true variability."))
    }
    stderrors <- apply(bootstrap_coefs, 2, sd, na.rm=TRUE)
    ci <- t(apply(bootstrap_coefs, 2, quantile, probs=c(a, 1-a), na.rm=TRUE))
    colnames(ci) <- c("lower", "upper")
    result <- list(est=est, se=stderrors, ci=ci, level=level, type="bootstrap", B=B)
  }
  return(result)
}
