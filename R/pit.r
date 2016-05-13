pit <- function(...) UseMethod("pit")

pit.default <- function(response, pred, distr=c("poisson", "nbinom"), distrcoefs, bins=10, ...){
  n <- length(pred)
  u <- seq(0, 1, length=bins+1)
  pit <- numeric(length(u))
  for(t in 1:n){
    P_x <- pdistr(response[t], meanvalue=pred[t], distr=distr, distrcoefs=distrcoefs)
    if(response[t]!=0){
      P_x_1 <- pdistr(response[t]-1, meanvalue=pred[t], distr=distr, distrcoefs=distrcoefs)
    }else{
      P_x_1 <- 0
    }
    pit <- pit + punif(u, P_x_1, P_x)/n
  }
  histo <- list(breaks=u, counts=diff(pit)*n, density=diff(pit)*bins, mids=(u[-(bins+1)]+u[-1])/2, xname="PIT", equidits=TRUE)
  class(histo) <- "histogram"
  plot_args <- modifyList(list(main="Non-randomized PIT histogram", xlab="Probability integral transform", ylab="Density", freq=FALSE, ylim=range(0, histo$density)), list(...)) #the default arguments can be overriden by those provided in the ... argument
  do.call("plot", args=c(list(x=histo), plot_args))
  #simconfint <- if(ci>0 && ci<1) (n/bins+c(-1,+1)*qnorm(1-(1-ci)/bins/2)*sqrt(n*(1/bins)*(1-1/bins)))/(n/bins) else NULL #simultaneous confidence band of level ci (normal approximation) for the histogram bars under the assumption of iid U(0,1) PIT values 
  #if(ci>0 && ci<1) abline(h=simconfint, lty="dashed", col=ci.col)
  abline(h=1, lty="dashed", col="blue")
}

pit.tsglm <- function(object, bins=10, ...){
  pit.default(response=object$response, pred=fitted(object), distr=object$distr, distrcoefs=object$distrcoefs, bins=bins, ...)
}
