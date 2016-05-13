tsglm.sim <- function(n, param=list(intercept=1, past_obs=NULL, past_mean=NULL, xreg=NULL), model=list(past_obs=NULL, past_mean=NULL, external=FALSE), xreg=NULL, link=c("identity", "log"), distr=c("poisson", "nbinom"), distrcoefs, fit, n_start=50){

  #Use information of a previous model fit if provided:
  if(!missing(fit)){
    tsglm.check(fit)
    #Some of the given arguments override the corresponding elements of the provided model fit in argument 'fit':    
    if(missing(n)){
      n <- fit$n_obs
    }
    if(missing(xreg)){
      xreg <- fit$xreg
    }
    if(missing(distr) || missing(distrcoefs)){
      distr <- fit$distr
      distrcoefs <- fit$distrcoefs
    }
    #Other arguments are preferably taken from the provided model fit in argument 'fit' to ensure a consistent model. If this is the case, it will be indicated to the user by a warning:
    if(!missing(param)) warning("Parameters provided in argument 'fit$coefficients' is used, the parameters in argument 'param' are ignored")
    paramvec <- fit$coefficients
    if(!missing(model)) warning("Model specification provided in argument 'fit$model' is used, the model specification in argument 'model' is ignored")
    model <- fit$model
    if(!missing(link)) warning("Link function provided in argument 'fit$link' is used, the link function in argument 'link' is ignored")
    link <- fit$link
  }else{
    if(missing(n)) stop("Argument 'n' is missing. Either provide this argument or an argument 'fit' with the result of a previous model fit")
  }
  
  #Check and modify arguments:
  link <- match.arg(link)
  distr <- match.arg(distr)
  checkdistr(distr=distr, distrcoefs=distrcoefs)
  model_names <- c("past_obs", "past_mean", "external")
  stopifnot(
    n%%1==0 & n>=0,
    n_start%%1==0 & n_start>=0,
    names(model) %in% model_names,
    model$past_obs%%1==0,
    model$past_mean%%1==0
  )
  model <- model[model_names]
  names(model) <- model_names
  if(is.null(xreg) || (is.matrix(xreg) && ncol(xreg)==0)) xreg <- matrix(0, nrow=n, ncol=0) else xreg <- as.matrix(xreg)
  p <- length(model$past_obs)
  P <- seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  p_max <- max(model$past_obs, 0)
  P_max <- seq(along=numeric(p_max)) #sequence 1:p_max if p_max>0 and NULL otherwise
  q <- length(model$past_mean)
  Q <- seq(along=numeric(q)) #sequence 1:q if q>0 and NULL otherwise
  q_max <- max(model$past_mean, 0)
  Q_max <- seq(along=numeric(q_max)) #sequence 1:q_max if q_max>0 and NULL otherwise
  r <- max(ncol(xreg), 0)
  R <- seq(along=numeric(r)) #sequence 1:r if r>0 and NULL otherwise
  if(p==0 & q>0) warning("Without dependence on past observations the dependence on past values of the linear predictor has no effect. Choose the model wisely.")
  stopifnot(
    nrow(xreg)==n,
    length(model$external)%in%c(r, 1, 0)
  )
  model$external <- as.logical(model$external)
  if(any(is.na(model$external))) stop("Argument 'model$external' could not be coerced to logical")
  if(length(model$external)==0) model$external <- rep(FALSE, r) #by default covariate has internal effect
  if(length(model$external)==1 && r!=1) rep(model$external, r) #if only one value is given this is used for all covariates
  if(!missing(fit)){ #use information of previous model fit if provided:
    param <- list( #transform parameter vector to a list
      intercept=paramvec[1],
      past_obs=paramvec[1+P],
      past_mean=paramvec[1+p+Q],
      xreg=paramvec[1+p+q+R]
    )
  }else{ #Check argument param for validity and modify it whenever necessary:
    param_names <- c("intercept", "past_obs", "past_mean", "xreg")
    stopifnot(
      names(param) %in% param_names,
      "intercept" %in% names(param),
      length(param$intercept)==1,
      length(param$past_obs)==p,
      length(param$past_mean)==q,
      length(param$xreg)==r
    )
    param <- param[param_names]
    names(param) <- param_names
    if(p==0) param$past_obs <- numeric() #allows to apply mathematical functions
    if(q==0) param$past_mean <- numeric() #allows to apply mathematical functions
    if(r==0) param$xreg <- numeric() #allows to apply mathematical functions
    if(is.null(model$past_obs)) model$past_obs <- seq(along=param$past_obs) #by default the parameters are for regression on the first few past observations
    if(is.null(model$past_mean)) model$past_mean <- seq(along=param$past_mean) #by default the parameters are for regression on the first few past means
  }
  
  #Check parameter restrictions:
  denom <- (1-sum(param$past_obs)-sum(param$past_mean))[[1]]    
  nu_stationary <- (param$intercept/denom)[[1]]
  tsglm.parametercheck(param, link=link)
  if(link=="identity") if(nu_stationary > 1e+09) stop("Too large mean to simulate from Poisson distribution, sum of parameters for regression on past observations and on past conditional means might be too close to one and/or intercept is too large")
     
  n_total <- n_start + n #total number of observations to be simulated (including burn-in)

  #Initialisation:
  if(n_start==0 & !missing(fit)){ #If simulation is based on a given fit and the length of the burn-in period is chosen to be 0, the the simulated observations are a direct continuation of the available observations.
    X_init <- fit$xreg[fit$n_obs-rev(Q_max)+1, , drop=FALSE]
    nu_init <- fit$linear.predictors[fit$n_obs-rev(Q_max)+1]
    y_init <- fit$ts[fit$n_obs-rev(P_max)+1]
  }else{
    X_init <- matrix(0, nrow=q_max+n_start, ncol=r) #the covariates during the burn-in period are set to zero because no other values are available
    nu_init <- rep(nu_stationary, q_max)
    y_init <- rdistr(n=p_max, meanvalue=g_inv(nu_stationary, link=link), distr=distr, distrcoefs=distrcoefs)
  }
  X <- rbind(X_init, xreg)    
  nu <- c(nu_init, numeric(n_total))  
  y <- c(y_init, integer(n_total))
    
  #Recursion:  
  for(t in 1:n_total){
    nu[t+q_max] <- param$intercept + sum(param$past_obs*trafo(y[(t-model$past_obs)+p_max], link=link)) + sum(param$past_mean*nu[(t-model$past_mean)+q_max]) + if(r>0){sum(param$xreg*X[t+q_max,]) - if(q>0){sum(param$past_mean*colSums(model$external*param$xreg*t(X[(t-model$past_mean)+q_max, , drop=FALSE])))}else{0}}else{0}
    y[t+p_max] <- rdistr(n=1, meanvalue=g_inv(nu[t+q_max], link=link), distr=distr, distrcoefs=distrcoefs)
  }    

  #Remove initialisation:
  y <- y[p_max+n_start+(1:n)]
  nu <- nu[q_max+n_start+(1:n)]
    
  y <- as.ts(y)
  nu <- as.ts(nu)
  xreg_effect <- as.ts(colSums(param$xreg*t(xreg)))
  result <- list(ts=y, linear.predictors=nu, xreg.effects=xreg_effect)
  return(result)
}
