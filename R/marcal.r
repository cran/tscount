marcal <- function(...) UseMethod("marcal")

marcal.default <- function(response, pred, distr=c("poisson", "nbinom"), distrcoefs, plot=TRUE, ...){
  xvalues <- min(response):max(response) #range of values could be extended in the future
  p_bar <- g_hat <- numeric(length(xvalues))
  n <- length(pred) 
  for(i in seq(along=xvalues)){
    for(t in 1:n){
      p_bar[i] <- p_bar[i] + pdistr(xvalues[i], meanvalue=pred[t], distr=distr, distrcoefs=distrcoefs)/n
    }
    g_hat[i] <- mean(response<=xvalues[i])
  }
  result <- list(x=xvalues, y=p_bar-g_hat)
  if(plot){
    plot_args <- modifyList(list(main="Marginal calibration plot", xlab="Threshold value", ylab="Diff. of pred. and emp. c.d.f"), list(...)) #the default arguments can be overriden by those provided in the ... argument
    do.call("plot", args=c(list(result$x, result$y, type="l"), plot_args))
    abline(h=0, lty=3)
    invisible(result)
  }else{
    return(result)
  }
}

marcal.tsglm <- function(object, plot=TRUE, ...){
  marcal.default(response=object$response, pred=fitted(object), distr=object$distr, distrcoefs=object$distrcoefs, plot=plot, ...)
}
