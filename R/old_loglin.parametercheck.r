#loglin.parametercheck <- function(param){
#  stopifnot(     
#    abs(param$past_obs)<1,
#    abs(param$past_mean)<1
#  )
#  sum_param_past <- abs(sum(param$past_obs)+sum(param$past_mean))
#  if(sum_param_past>=1) stop(paste("Parameters are outside the stationary region, absolute sum of parameters for\nregression on past observations and on past conditional means is", sum_param_past, "> 1"))
#  return(TRUE)
#}
 