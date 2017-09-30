tsglm.loglik <- function(paramvec, model, ts, xreg, link, score=FALSE, info=c("none", "score", "hessian", "sandwich"), condmean=NULL, from=1, init.method=c("marginal", "iid", "firstobs", "zero"), init.drop=FALSE){
  #Conditional (quasi) log-likelihood function, score function and information matrix of a count time series following generalised linear models
  
  ##############                  
  #Checks and preparations:
  init.method <- match.arg(init.method)
  n <- length(ts)
  p <- length(model$past_obs)
  p_max <- max(model$past_obs, 0)
  q <- length(model$past_mean)
  q_max <- max(model$past_mean, 0)
  r <- max(ncol(xreg), 0)
  R <- seq(along=numeric(r)) #sequence 1:r if r>0 and NULL otherwise
  info <- match.arg(info)
  if(!score & info!="none"){
    score <- TRUE
    warning("Information matrix cannot be calculated without score vector. Argument 'score'\nis set to TRUE.")
  }
  derivatives <- if(!score) "none" else if(info %in% c("hessian", "sandwich")) "second" else "first"
  parameternames <- tsglm.parameternames(model=model, xreg=xreg)
  startfrom <- ifelse(init.drop, p_max+1, 1) #first time point which is considered for the final result
  n_effective <- ifelse(init.drop, n-p_max, n) #effective number of observations considered for the final result
  ##############
  
  condmean <- tsglm.condmean(paramvec=paramvec, model=model, ts=ts, xreg=xreg, link=link, derivatives=derivatives, condmean=condmean, from=from, init.method=init.method)
  #Load objects and remove initialisation if necessary:
  z <- condmean$z[p_max+(startfrom:n)]
  nu <- condmean$nu[q_max+(startfrom:n)]
  if(derivatives %in% c("first", "second")){
    partial_nu <- condmean$partial_nu[q_max+(startfrom:n), , drop=FALSE]
    partial_lambda <- g_inv_1st(nu, link=link)*partial_nu
  }
  if(derivatives == "second"){
    partial2_nu <- condmean$partial2_nu[q_max+(startfrom:n), , , drop=FALSE]
    outerprod_partial2_nu <- array(NA, dim=c(n_effective, 1+p+q+r, 1+p+q+r), dimnames=list(NULL, parameternames, parameternames))
    outerprod_partial2_nu[] <- if(p+q+r > 0) aperm(sapply(1:n_effective, function(i) partial_nu[i,]%*%t(partial_nu[i,]), simplify="array"), c(3,1,2)) else array(partial_nu[,1]^2, dim=c(n_effective,1,1))
    partial2_lambda <- g_inv_2nd(nu, link=link)*outerprod_partial2_nu + g_inv_1st(nu, link=link)*partial2_nu
  }
  lambda <- g_inv(nu, link=link)
  y <- trafo_inv(z, link=link)   
  loglik_t <- ifelse(lambda>0, y*log(lambda)-lambda, -Inf)
  loglik <- sum(loglik_t)
  scorevec <- NULL
  if(score){
    scorevec_t <- (y/lambda-1) * partial_lambda
    scorevec <- colSums(scorevec_t)
  }
  outerscoreprod <- NULL
  infomat <- NULL
  if(info != "none"){
    if(info %in% c("score", "sandwich")){
      outerscoreprod <- array(NA, dim=c(n_effective, 1+p+q+r, 1+p+q+r), dimnames=list(NULL, parameternames, parameternames))
      outerscoreprod[] <- if(p+q+r > 0) aperm(sapply(1:n_effective, function(i) partial_lambda[i,]%*%t(partial_lambda[i,]), simplify="array"), c(3,1,2)) else array(partial_lambda[,1]^2, dim=c(n_effective,1,1))
      infomat <- infomat_score <- apply(1/lambda*outerscoreprod, c(2,3), sum)
    }
    if(info %in% c("hessian", "sandwich")){
      hessian_t <- aperm((-y/lambda^2) * replicate(1+p+q+r, partial_lambda) * aperm(replicate(1+p+q+r, partial_lambda), perm=c(1,3,2)), perm=c(2,3,1)) + rep((y/lambda-1), each=(1+p+q+r)^2) * aperm(partial2_lambda, perm=c(2,3,1))
      infomat <- infomat_hessian <- -apply(hessian_t, c(1,2), sum)
    }
    if(info == "sandwich"){
      infomat <- infomat_hessian %*% invertinfo(infomat_score, stopOnError=TRUE)$vcov %*% infomat_hessian
      outerscoreprod <- NULL
    }
    dimnames(infomat) <- list(parameternames, parameternames) 
  }
  result <- list(loglik=loglik, score=scorevec, info=infomat, outerscoreprod=outerscoreprod, nu=nu)
  return(result)
}
