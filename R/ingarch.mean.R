ingarch.mean <- function(intercept, past_obs=NULL, past_mean=NULL){
#Theoretical marginal mean of an INGARCH(p,q) process
##############################
  tsglm.parametercheck(param=list(intercept=intercept, past_obs=past_obs, past_mean=past_mean, xreg=NULL), link="identity")    
  result <- (intercept/(1-sum(past_mean)-sum(past_obs)))[[1]]
  return(result)
}
