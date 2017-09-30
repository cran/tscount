scoring <- function(...) UseMethod("scoring")

scoring.default <- function(response, pred, distr=c("poisson", "nbinom"), distrcoefs, individual=FALSE, cutoff=1000, ...){
  n <- length(pred) 
  logarithmic <- quadratic <- spherical <- rankprob <- dawseb <- normsq <- sqerror <- numeric(n) #scores
  for(t in 1:n){
    #auxiliary objects, which are overwritten in each step:
      y <- response[t]
      mu <- pred[t]
      sigma <- sddistr(meanvalue=mu, distr=distr, distrcoefs=distrcoefs)
      p_y <- ddistr(y, meanvalue=mu, distr=distr, distrcoefs=distrcoefs)
      quadrat_p <- sum(ddistr(0:cutoff, meanvalue=mu, distr=distr, distrcoefs=distrcoefs)^2)
    #computation of the scores:
      logarithmic[t] <- - log(p_y)
      quadratic[t] <- - 2*p_y + quadrat_p
      spherical[t] <- - p_y/sqrt(quadrat_p)
      rankprob[t] <- sum((pdistr(0:cutoff, meanvalue=mu, distr=distr, distrcoefs=distrcoefs) - as.integer(y <= 0:cutoff))^2)
      sqerror[t] <- (y-mu)^2
      normsq[t] <- sqerror[t]/sigma^2 
      dawseb[t] <- normsq[t] + 2*log(sigma)
  }
  result <- data.frame(
    logarithmic=logarithmic,
    quadratic=quadratic,
    spherical=spherical,
    rankprob=rankprob,
    dawseb=dawseb,
    normsq=normsq,
    sqerror=sqerror
  )
  if(!individual) result <- apply(result, 2, mean)
  return(result)
}

scoring.tsglm <- function(object, individual=FALSE, cutoff=1000, ...){
  scoring.default(response=object$response, pred=fitted(object), distr=object$distr, distrcoefs=object$distrcoefs, individual=individual, cutoff=cutoff, ...)
}
