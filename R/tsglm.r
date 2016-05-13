tsglm <- function(ts, model=list(past_obs=NULL, past_mean=NULL, external=NULL), xreg=NULL, link=c("identity", "log"), distr=c("poisson", "nbinom"), ...){
  link <- match.arg(link)  
  distr <- match.arg(distr)
  cl <- match.call()
  #Estimating the mean structure:
  meanfit <- mean.fit(ts=ts, model=model, xreg=xreg, link=link, ...)
  if(length(meanfit$coefficients)==0){ #if no final estimation is done, then the function returns a list with less elements and without the class 'tsglm'
    result <- list(start=meanfit$start, call=cl, n_obs=meanfit$n_obs, ts=meanfit$ts, model=meanfit$model, xreg=meanfit$xreg, link=link)
    return(result)
  }
  #Estimating the distribution:
  disfit <- distr.fit(meanfit, distr=distr)
  info.matrix_corrected <- if(is.null(meanfit$outerscoreprod)) NULL else apply(as.numeric(1/meanfit$fitted.values + disfit$sigmasq)*meanfit$outerscoreprod, c(2,3), sum)
  loglik <- sum(ddistr(ts, meanvalue=meanfit$fitted.values, distr=disfit$distr, distrcoefs=disfit$distrcoefs, log=TRUE))
  result <- c(
    list(coefficients=meanfit$coefficients, start=meanfit$start, residuals=meanfit$residuals, fitted.values=meanfit$fitted.values, linear.predictors=meanfit$linear.predictors, response=meanfit$response, logLik=loglik, score=meanfit$score, info.matrix=meanfit$info.matrix, info.matrix_corrected=info.matrix_corrected, call=cl, n_obs=meanfit$n_obs, n_eff=meanfit$n_eff, ts=meanfit$ts, model=meanfit$model, xreg=meanfit$xreg, link=link), #an extract of the object meanfit with an additional information matrix corrected for the fitted conditional distribution
    disfit
  )
  class(result) <- c("tsglm")
  return(result)
}
