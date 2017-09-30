residuals.tsglm <- function(object, type=c("response", "pearson", "anscombe"), ...){
  type <- match.arg(type)
  result <- switch(type,
    "response"=object$residuals,  
    "pearson"=object$residuals/sddistr(meanvalue=fitted(object), distr=object$distr, distrcoefs=object$distrcoefs),
    "anscombe"=ardistr(response=object$response, meanvalue=fitted(object), distr=object$distr, distrcoefs=object$distrcoefs)
  )
  return(result)  
}
