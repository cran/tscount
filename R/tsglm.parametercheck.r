tsglm.parametercheck <- function(param, link=c("identity", "log"), stopOnError=TRUE, silent=TRUE){
  #Check parameter vector of a count time series following GLMs

  ##############
  #Checks and preparations:
  link <- match.arg(link)
  
  if(link == "identity") parametercheck <- function(param){
    stopifnot(     
      param$intercept>0,
      param$past_obs>=0,
      param$past_mean>=0,
      param$xreg>=0
    )
    sum_param_past <- sum(param$past_obs)+sum(param$past_mean)
    if(sum_param_past>=1) stop(paste("Parameters are outside the stationary region, sum of parameters for regression\non past observations and on past conditional means is", sum_param_past, "> 1"))
    return(TRUE)
  }  

  if(link == "log") parametercheck <- function(param){
    stopifnot(     
      abs(param$past_obs)<1,
      abs(param$past_mean)<1
    )
    sum_param_past <- abs(sum(param$past_obs)+sum(param$past_mean))
    if(sum_param_past>=1) stop(paste("Parameters are outside the stationary region, absolute sum of parameters for\nregression on past observations and on past conditional means is", sum_param_past, "> 1"))
    return(TRUE)
  }
  
  if(stopOnError){  
    result <- parametercheck(param)  
  }else{
    result <- try(parametercheck(param), silent=silent)
    if(class(result)=="try-error") result <- FALSE
  } 
  return(result)
}
