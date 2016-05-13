tsglm.condmean <- function(paramvec, model, ts, xreg, link, derivatives=c("none", "first", "second"), condmean=NULL, from=1, init.method=c("marginal", "iid", "firstobs", "zero")){
  #Recursion for the linear predictor (which is the conditional mean for the identity link) and its derivatives of a count time series following generalised linear models
 
  ##############
  #Checks and preparations:
  n <- length(ts)
  p <- length(model$past_obs)
  P <- seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  p_max <- max(model$past_obs, 0)
  q <- length(model$past_mean)
  Q <- seq(along=numeric(q)) #sequence 1:q if q>0 and NULL otherwise
  q_max <- max(model$past_mean, 0)
  Q_max <- seq(along=numeric(q_max))
  r <- max(ncol(xreg), 0)
  R <- seq(along=numeric(r)) #sequence 1:r if r>0 and NULL otherwise
  parameternames <- tsglm.parameternames(model=model, xreg=xreg)
  derivatives <- match.arg(derivatives)
  init.method <- match.arg(init.method)
  param <- list( #transform parameter vector to a list
    intercept=paramvec[1],
    past_obs=paramvec[1+P],
    past_mean=paramvec[1+p+Q],
    xreg=paramvec[1+p+q+R]
  )    
  if(!is.null(condmean)){ #If the output of a previous call is provided, the recursion starts from t=from. Else initialisation of all objects is necessary and the recursion starts from t=1.
    times <- if(from <= n) from:n else NULL
    #Load objects:
    z <- condmean$z
    nu <- condmean$nu
    if(derivatives %in% c("first", "second")) partial_nu <- condmean$partial_nu
    if(derivatives == "second") partial2_nu <- condmean$partial2_nu
#########include checks if argument condmean is sufficient for further calculations
  }else{  
    times <- 1:n

  # # # # # # #
  #Initialisation:    
    if(init.method == "marginal"){ #initialisation by stationary solution (and its partial derivatives)
      denom <- (1-sum(param$past_obs)-sum(param$past_mean))[[1]]    
      nu_stationary <- (param$intercept/denom)[[1]]
      nu <- c(rep(nu_stationary, q_max), numeric(n))  
      z <- c(rep(nu_stationary, p_max), trafo(ts, link=link))
      if(derivatives %in% c("first", "second")){
        #Vector of first partial derivatives of nu with respect to the parameters:
        partial_nu <- matrix(0, nrow=n+q_max, ncol=1+p+q+r)
        partial_nu[Q_max, 1] <- 1/denom #intercept
        partial_nu[Q_max, 1+P] <- param$intercept/denom^2 #past_obs
        partial_nu[Q_max, 1+p+Q] <- param$intercept/denom^2 #past_mean
        #derivatives with respect to the regressor coefficients are zero (which is the default)
        if(derivatives == "second"){
          #Matrix of second partial derivatives of nu with respect to the parameters:
          partial2_nu <- array(0, dim=c(n+q_max, 1+p+q+r, 1+p+q+r))  
          partial2_nu[Q_max, 1, 1+c(P,p+Q)] <- partial2_nu[Q_max, 1+c(P,p+Q), 1] <- 1/denom^2
          partial2_nu[Q_max, 1+c(P,p+Q), 1+c(P,p+Q)] <- 2*param$intercept/denom^3
          #derivatives with respect to the regressor coefficients are zero
        }
      }
    }
    if(init.method == "iid"){ #initialisation under iid assumption:
      nu <- c(rep(param$intercept, q_max), numeric(n))  
      z <- c(rep(param$intercept, p_max), trafo(ts, link=link))
      if(derivatives %in% c("first", "second")){
        partial_nu <- matrix(0, nrow=n+q_max, ncol=1+p+q+r)
        partial_nu[Q_max, 1] <- 1 #intercept
        if(derivatives == "second"){
          partial2_nu <- array(0, dim=c(n+q_max, 1+p+q+r, 1+p+q+r))  
        }
      }
    }
    if(init.method == "firstobs"){ #initialisation by the first observation:
      firstobs <- ts[1]
      nu <- c(rep(trafo(firstobs, link=link), q_max), numeric(n))  
      z <- c(rep(trafo(firstobs, link=link), p_max), trafo(ts, link=link))
      if(derivatives %in% c("first", "second")){
        partial_nu <- matrix(0, nrow=n+q_max, ncol=1+p+q+r)
        if(derivatives == "second"){
          partial2_nu <- array(0, dim=c(n+q_max, 1+p+q+r, 1+p+q+r))  
        }
      }      
    }
    if(init.method == "zero"){ #initialisation by value zero:
      nu <- c(rep(0, q_max), numeric(n))  
      z <- c(rep(trafo(0, link=link), p_max), trafo(ts, link=link))
      if(derivatives %in% c("first", "second")){
        partial_nu <- matrix(0, nrow=n+q_max, ncol=1+p+q+r)
        if(derivatives == "second"){
          partial2_nu <- array(0, dim=c(n+q_max, 1+p+q+r, 1+p+q+r))  
        }
      }
    }    
  }
  X <- matrix(0, nrow=q_max+n, ncol=r) #regressors are initalised by zero
  X[q_max+(1:n), ] <- xreg
  # # # # # # #
  ##############
  
  ##############
  #Recursion:
  
  # # # # # # #
  #Conditional mean:
  for(t in times){
    nu[t+q_max] <- param$intercept + sum(param$past_obs*z[(t-model$past_obs)+p_max]) + sum(param$past_mean*nu[(t-model$past_mean)+q_max]) + if(r>0){sum(param$xreg*X[t+q_max, ]) - if(q>0){sum(param$past_mean*colSums(model$external*param$xreg*t(X[(t-model$past_mean)+q_max, , drop=FALSE])))}else{0}}else{0}  
  }
  result <- list(z=z, nu=nu)
  # # # # # # #
  
  # # # # # # #
  #First derivatives:    
  if(derivatives %in% c("first", "second")){
    for(t in times){
      partial_nu[t+q_max, 1] <- 1 + sum(param$past_mean*partial_nu[(t-model$past_mean)+q_max, 1]) #intercept
      if(p>0) partial_nu[t+q_max, 1+P] <- z[t-model$past_obs+p_max] + (if(q>0){t(param$past_mean) %*% partial_nu[(t-model$past_mean)+q_max, 1+P, drop=FALSE]}else{numeric(p)}) #past_obs    
      if(q>0) partial_nu[t+q_max, 1+p+Q] <- nu[t-model$past_mean+q_max] + t(param$past_mean) %*% partial_nu[(t-model$past_mean)+q_max, 1+p+Q, drop=FALSE] - (if(r>0){param$past_mean*colSums(model$external*param$xreg*t(X[(t-model$past_mean)+q_max, , drop=FALSE]))}else{numeric(q)}) #past_mean
      if(r>0) partial_nu[t+q_max, 1+p+q+R] <- (if(q>0){colSums(param$past_mean*partial_nu[(t-model$past_mean)+q_max, 1+p+q+R, drop=FALSE]) - model$external*colSums(param$past_mean*X[(t-model$past_mean)+q_max, , drop=FALSE])}else{numeric(r)}) + X[t+q_max, ] #covariates
    }
    dimnames(partial_nu)[[2]] <- if(p==0 & q==0 & r==0) list(parameternames) else parameternames
    result <- c(result, list(partial_nu=partial_nu))
  # # # # # # #
  
  # # # # # # #
  #Second derivative:
    if(derivatives == "second"){
      for(t in times){
        partial2_nu[t+q_max, , ] <- apply(param$past_mean*partial2_nu[t+q_max-model$past_mean, , , drop=FALSE], c(2,3), sum)
        partial2_nu[t+q_max, 1+p+Q, 1+p+Q] <- partial2_nu[t+q_max, 1+p+Q, 1+p+Q] + (partial_nu[t+q_max-model$past_mean, 1+p+Q] + t(partial_nu[t+q_max-model$past_mean, 1+p+Q]))/2 #from the formula we would only need the first part of the last summand, but in this case our matrix is not symmetrical, so we add this average
        partial2_nu[t+q_max, 1+p+Q, 1+p+q+R] <- partial2_nu[t+q_max, 1+p+Q, 1+p+q+R] + partial_nu[t+q_max-model$past_mean, 1+p+q+R] - X[t+q_max-model$past_mean,]
        partial2_nu[t+q_max, 1+p+q+R, 1+p+Q] <- t(partial2_nu[t+q_max, 1+p+Q, 1+p+q+R])
      }
      dimnames(partial2_nu)[[2]] <- dimnames(partial2_nu)[[3]] <- if(p==0 & q==0 & r==0) list(parameternames) else parameternames
      result <- c(result, list(partial2_nu=partial2_nu))  
    }    
  }
  # # # # # # #  
  ##############
    
  return(result)
} 
