QIC <- function(...) UseMethod("QIC")

QIC.tsglm <- function(object, ...){
  qloglik <- sum(ddistr(object$ts, meanvalue=fitted(object), distr="poisson", log=TRUE)) #Q
  V <- vcov(object)
  O <- tsglm.loglik(paramvec=coef(object), model=object$model, ts=object$ts, xreg=object$xreg, link=object$link, score=TRUE, info="hessian", ...)$info #it would be more consistent to choose arguments 'init.method' and 'init.drop' like in the original function call
  dof <- sum(diag(O%*%V))
  result <- -2*qloglik+2*dof
  return(result)
}
