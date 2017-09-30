#The R code in this file produces the RData files which are needed to compile the Rnw file of the article. The RData files are part of the submission and do not need to be generated again. Running this R code is only feasible on a high performance computing cluster.

#submitR parallel tscount-computations.r walltime=48:00:00 memory=1028 ntasks=51 email=liboschik

#screen
#interactive walltime=96:00:00 memory=1024 ncpus=11

#ntasks <- 10

#setwd("Forschung/Zaehldatenzeitreihen/Package/Journal of Statistical Software/Results")

#Load package:
library(tscount)

#Setting up parallel computation:
library(Rmpi)
library(parallel)
ntasks <- detectCores() - 1
clust <- makePSOCKcluster(ntasks, type="MPI") #for use on local Windows PC
#clust <- makeCluster(ntasks, type="MPI")
setDefaultCluster(cl=clust)
print(clust)



# # # # # # #
# Seatbelts #
# # # # # # #

timeseries <- Seatbelts[, "VanKilled"]
regressors <- cbind(PetrolPrice=Seatbelts[, c("PetrolPrice")],
                    linearTrend=seq(along=timeseries)/12)
timeseries_until1981 <- window(timeseries, end=1981+11/12)
regressors_until1981 <- window(regressors, end=1981+11/12)
set.seed(1622)
seatbeltsfit <- tsglm(ts=timeseries_until1981,
  model=list(past_obs=c(1, 12)), link="log", distr="pois",
  xreg=regressors_until1981)
durationse <- system.time(seatbeltssummary <- summary(seatbeltsfit, B=500, parallel=TRUE))
warningse <- warnings()
save(seatbeltssummary, durationse, warningse, file="seatbelts.RData")



# # # # # # # # #
# Campylobacter #
# # # # # # # # #

interventions <- interv_covariate(n=length(campy), tau=c(84, 100), delta=c(1, 0))
campymodel <- list(past_obs=1, past_mean=13)
set.seed(1622)
campyfit_nbin <- tsglm(campy, model=campymodel, xreg=interventions, dist="nbinom")
durationse <- system.time(campyse <- se(campyfit_nbin, B=500))
warningse <- warnings()
save(campyse, durationse, warningse, file="campy.RData")


# # # # # # # #
# Covariates  #
# # # # # # # #

#How to choose the intercept parameter for the log-linear model such that the marginal mean is about 4? Confirm try-and-error solution by a simulation:
#mean(tsglm.sim(n=1000000, param=list(intercept=0.65, past_obs=0.3, past_mean=0.2), model=list(past_obs=1, past_mean=1), link="log")$ts)

covariate_sim1 <- function(seed, covariate_type, covariate_factor=1, link="log", n=100, param=list(intercept=ifelse(link=="identity", 2, 0.65), past_obs=0.3, past_mean=0.2), plot=FALSE, fitonly=FALSE, ...){
  require(tscount)
  if(!missing(seed)) set.seed(seed)
  linear <- (1:n)/n
  if(covariate_type=="linear") covariate <- linear
  if(covariate_type=="quadratic") covariate <- linear^2
  if(covariate_type=="sine_fixednumber") covariate <- (sin(2*pi*linear*5)+1)/2 #five full periods
  if(covariate_type=="sine_fixedlength") covariate <- (sin(2*pi*linear*(n/20))+1)/2 #period of length 20 time units
  if(covariate_type=="interv_so") covariate <- interv_covariate(n, tau=floor(n/2), delta=0)
  if(covariate_type=="interv_ts") covariate <- interv_covariate(n, tau=floor(n/2), delta=0.8) 
  if(covariate_type=="interv_ls") covariate <- interv_covariate(n, tau=floor(n/2), delta=1)  
  if(covariate_type=="garch"){
    require(fGarch)
    covariate <- as.numeric(garchSim(spec=garchSpec(model=list(omega=0.002, alpha=0.1, beta=0.8, mu=0.5), cond.dist="norm"), n=n))
  }
  if(covariate_type=="iid_pois") covariate <- rpois(n, lambda=0.5)
  if(covariate_type=="iid_exp") covariate <- rexp(n, rate=2)
  if(covariate_type=="iid_norm") covariate <- rnorm(n, mean=0.5, sd=0.2)
  if(covariate_type=="iid_chisq") covariate <- rchisq(n, df=0.5)
  covariate <- pmax(covariate, 0) #ensures that covariate is non-negative
  model <- list(past_obs=1, past_mean=1, external=TRUE)
  param <- c(param, list(xreg=covariate_factor*param$intercept))
  timser <- tsglm.sim(n, param=param, model=model, xreg=covariate, link=link)$ts
  fit <- try(tsglm.meanfit(ts=timser, model=model, xreg=covariate, link=link, ...))
  if(class(fit)=="try-error") return(rep(NA, 11))
  if(fitonly) return(fit)
  if(plot){
    plot(timser, type="o", ylim=c(0, max(timser)))
    if(link=="identity") lines(covariate*param$xreg, lty="dashed")
    if(link=="log") lines(exp(covariate*param$xreg), lty="dashed")
    print(summary(fit))
  }
  result <- c(final=fit$coefficients, start=fit$start, loglik=fit$logLik, optimtotal=fit$final$counts[1], optimouter=fit$final$outer.iterations)
  return(result)
}

covariate_sim <- function(N, settings, covariate_types=c("linear", "quadratic", "sine_fixednumber", "sine_fixedlength", "interv_so", "interv_ts", "interv_ls", "garch", "iid_pois", "iid_exp", "iid_norm", "iid_chisq"), seed=1027, ...){
  if(!is.null(seed)) set.seed(seed)
  seeds <- sample(1e+8, size=N)
  durations_covariate <- numeric(length(covariate_types))
  names(durations_covariate) <- covariate_types
  estimates_covariate <- vector("list", length(covariate_types))
  names(estimates_covariate) <- covariate_types
  for(covariate_type in covariate_types){
    durations_covariate[[covariate_type]] <- system.time(estimates_covariate[[covariate_type]] <- parSapply(cl=clust, seeds, covariate_sim1, covariate_type=covariate_type, covariate_factor=settings$covariate_factor, n=settings$n, link=settings$link, ...))[3]
  }
  result <- list(
    settings=settings,
    estimates=estimates_covariate,
    durations=durations_covariate
  )
  print(c(duration=sum(durations_covariate)))
  return(result)
}

##Note that there is one erroneous simulation run for the chi^2 covariate. Consider to make the simulation more robust before running it another time. 

N <- 200 #number of simulated time series

covariate_n100_id <- covariate_sim(N=N, settings=list(n=100, link="identity", covariate_factor=2))
covariate_n100_log <- covariate_sim(N=N, settings=list(n=100, link="log", covariate_factor=1.5))
covariate_n500_id <- covariate_sim(N=N, settings=list(n=500, link="identity", covariate_factor=2))
covariate_n500_log <- covariate_sim(N=N, settings=list(n=500, link="log", covariate_factor=1.5))
covariate_n1000_id <- covariate_sim(N=N, settings=list(n=1000, link="identity", covariate_factor=2))
covariate_n1000_log <- covariate_sim(N=N, settings=list(n=1000, link="log", covariate_factor=1.5))
covariate_n2000_id <- covariate_sim(N=N, settings=list(n=2000, link="identity", covariate_factor=2))
covariate_n2000_log <- covariate_sim(N=N, settings=list(n=2000, link="log", covariate_factor=1.5))

save(
  covariate_n100_id, covariate_n100_log,
  covariate_n500_id, covariate_n500_log,
  covariate_n1000_id, covariate_n1000_log,
  covariate_n2000_id, covariate_n2000_log,
file="covariates.RData")



#####Choice of start estimation:
###Does the estimation depend on the start estimation? Does the estimation become better with the true parameters as starting values?
#
#covariate_GLMstart <- covariate_sim(N=N, settings=list(n=1000, covariate_factor=2, link="identity"), covariate_types=c("linear", "sine_fixednumber", "interv_ls", "garch"), seed=1027, start.control=list(method="GLM"))
#covariate_ARMAstart <- covariate_sim(N=N, settings=list(n=1000, covariate_factor=2, link="identity"), covariate_types=c("linear", "sine_fixednumber", "interv_ls", "garch"), seed=1027, start.control=list(method="CSS"))
#covariate_truestart <- covariate_sim(N=N, settings=list(n=1000, covariate_factor=2, link="identity"), covariate_types=c("linear", "sine_fixednumber", "interv_ls", "garch"), seed=1027, start.control=list(method="fixed", intercept=0.5*4, past_obs=0.3, past_mean=0.2, xreg=0.5*4*2))
#covariate_fixedstart <- covariate_sim(N=N, settings=list(n=1000, covariate_factor=2, link="identity"), covariate_types=c("linear", "sine_fixednumber", "interv_ls", "garch"), seed=1027, start.control=list(method="fixed"))
#save(
#  covariate_GLMstart,
#  covariate_ARMAstart,
#  covariate_truestart,
#  covariate_fixedstart,
#file="covariates_start.RData")
#
#load("covariates_start.RData")
#covariate_type <- "linear" #or any of the other considered types of covariates
#estimatesARMA <- covariate_ARMAstart$estimates[[covariate_type]]
#estimatesGLM <- covariate_GLMstart$estimates[[covariate_type]]
#estimatestrue <- covariate_truestart$estimates[[covariate_type]]
#estimatesfixed <- covariate_fixedstart$estimates[[covariate_type]]
#
#par(mfrow=c(3,2))
#yrange <- range(estimatesARMA[c(4,8),], estimatesGLM[c(4,8),], estimatestrue[c(4,8),], estimatesfixed[c(4,8),]) 
#highlight <- estimatesARMA[4, ] > 5.5
#pchs <- ifelse(highlight, 4, 1)
#cols <- ifelse(highlight, "red", "black")
#plot(colSums(estimatesARMA[c(2:3)+4, ]), estimatesARMA[4+4, ], xlab="beta_1 + alpha_1", ylab="eta_1", xlim=c(0.2,0.7), ylim=yrange, pch=pchs, col=cols, main="ARMA")
#abline(h=4); abline(v=0.5)
#plot(colSums(estimatesARMA[c(2:3), ]), estimatesARMA[4, ], xlab="beta_1 + alpha_1", ylab="eta_1", xlim=c(0.2,0.7), ylim=yrange, pch=pchs, col=cols, main="ARMAstart")
#abline(h=4); abline(v=0.5)
#plot(colSums(estimatesGLM[c(2:3)+4, ]), estimatesGLM[4+4, ], xlab="beta_1 + alpha_1", ylab="eta_1", xlim=c(0.2,0.7), ylim=yrange, pch=pchs, col=cols, main="GLM")
#abline(h=4); abline(v=0.5)
#plot(colSums(estimatesGLM[c(2:3), ]), estimatesGLM[4, ], xlab="beta_1 + alpha_1", ylab="eta_1", xlim=c(0.2,0.7), ylim=yrange, pch=pchs, col=cols, main="GLMstart")
#abline(h=4); abline(v=0.5)
#plot(colSums(estimatestrue[c(2:3), ]), estimatestrue[4, ], xlab="beta_1 + alpha_1", ylab="eta_1", xlim=c(0.2,0.7), ylim=yrange, pch=pchs, col=cols, main="truestart")
#abline(h=4); abline(v=0.5)
#plot(colSums(estimatesfixed[c(2:3), ]), estimatesfixed[4, ], xlab="beta_1 + alpha_1", ylab="eta_1", xlim=c(0.2,0.7), ylim=yrange, pch=pchs, col=cols, main="fixedstart")
#abline(h=4); abline(v=0.5)
#
#par(mfrow=c(1,2))
#plot(estimatesGLM[9, ], estimatesARMA[9, ], col=cols, pch=pchs, xlab="Loglik GLMstart", ylab="Loglik ARMAstart")
#plot(estimatestrue[9, ], estimatesfixed[9, ], col=cols, pch=pchs, xlab="Loglik truestart", ylab="Loglik fixedstart")
#
#
###How quick is the final maximum likelihood estimation for the different starting values?
#par(mfrow=c(2,2))
#hist(estimatesARMA[10, ], xlim=c(0,1200), freq=FALSE, main="ARMAstart")
#hist(estimatesGLM[10, ], xlim=c(0,1200), freq=FALSE, main="GLMstart")
#hist(estimatestrue[10, ], xlim=c(0,1200), freq=FALSE, main="truestart")
#hist(estimatesfixed[10, ], xlim=c(0,1200), freq=FALSE, main="fixedstart")
#sum(estimatesARMA[10, ])
#sum(estimatesGLM[10, ])
#sum(estimatestrue[10, ])
#sum(estimatesfixed[10, ])
##ARMA leads to a quicker optimisation than GLM, but has the problem of not converging every time discussed above. GLM does not have this problem, but takes even longer than using arbitary  starting values.
#
#
###How does the likelihood look like for a bad and for a nice situation?
#which(highlight)
#set.seed(1027)
#seeds <- sample(1e+8, size=200)
#summary(t(estimatesARMA[1:4,])) #choose range for evaluating the likelihood
##Bad situation:
#system.time(badfit <- covariate_sim1(seed=seeds[5], covariate_type="linear", covariate_factor=2, link="identity", n=1000, param=list(intercept=0.5*4, past_obs=0.3, past_mean=0.2), fitonly=TRUE))
#coef(badfit)
#system.time(badfit2 <- covariate_sim1(seed=seeds[5], covariate_type="linear", covariate_factor=2, link="identity", n=1000, param=list(intercept=0.5*4, past_obs=0.3, past_mean=0.2), fitonly=TRUE, final.control=list(optim.method="Nelder-Mead")))
#coef(badfit2)
#interceptvals <- coef(badfit)[1] #seq(from=1, to=3.5, length.out=10)
#past_obsvals <- coef(badfit)[2] #seq(from=0.1, to=0.4, length.out=40)
#past_meanvals <- seq(from=0, to=0.7, length.out=30)
#xregvals <- seq(from=0, to=8, length.out=30)
#pargrid <- expand.grid(
#  intercept=interceptvals,
#  past_obs=past_obsvals,
#  past_mean=past_meanvals,
#  xreg=xregvals
#)
#system.time(badlogliks <- apply(pargrid, 1, function(x) tscount:::ingarch.loglik(paramvec=x, model=badfit$model, xreg=badfit$xreg, ts=badfit$ts)$loglik))
##Nice situation:
#nicefit <- covariate_sim1(seed=seeds[6], covariate_type="linear", covariate_factor=2, link="identity", n=1000, param=list(intercept=0.5*4, past_obs=0.3, past_mean=0.2), fitonly=TRUE)
#coef(nicefit)
#interceptvals <- coef(nicefit)[1] #seq(from=1, to=3.5, length.out=10)
#past_obsvals <- coef(nicefit)[2] #seq(from=0.1, to=0.4, length.out=40)
#past_meanvals <- seq(from=0, to=0.7, length.out=30)
#xregvals <- seq(from=0, to=8, length.out=30)
#pargrid <- expand.grid(
#  intercept=interceptvals,
#  past_obs=past_obsvals,
#  past_mean=past_meanvals,
#  xreg=xregvals
#)
#system.time(nicelogliks <- apply(pargrid, 1, function(x) tscount:::ingarch.loglik(paramvec=x, model=nicefit$model, xreg=nicefit$xreg, ts=nicefit$ts)$loglik))
##Plot:
#par(mfrow=c(1,2))
#contour(past_meanvals, xregvals, matrix(badlogliks, ncol=length(xregvals)), levels=c(seq(5000, 7130, by=10)), main="Bad")
#points(coef(badfit)[3], coef(badfit)[4], pch=4)
#points(badfit$start[3], badfit$start[4], pch=1)
#contour(past_meanvals, xregvals, matrix(nicelogliks, ncol=length(xregvals)), levels=c(seq(5000, 7130, by=10)), main="Nice")
#points(coef(nicefit)[3], coef(nicefit)[4], pch=4)
#points(nicefit$start[3], nicefit$start[4], pch=1)
##The likelihood does not look very different in the two sitautions. In the bad situation the BFGS algorithm does not find the global optimum, even for a lower value for reltol. Using the Nelder-Mead algorithm helps but takes more than four times the computation time. Just using the fixed starting values is sufficient for solving the problem.



# # # # # # # # # # # # #
# Dispersion parameter  #
# # # # # # # # # # # # #

distrcoef_sim1 <- function(seed, size, link, n, param=list(intercept=ifelse(link=="identity", 0.5*4, 0.5*log(4)), past_obs=0.3, past_mean=0.2), ...){
  require(tscount)
  if(!missing(seed)) set.seed(seed)
  model <- list(past_obs=1, past_mean=1)
  if(size==Inf){
    timser <- tsglm.sim(n=n, param=param, model=model, link=link, distr="poisson")$ts
  }else{
    timser <- tsglm.sim(n=n, param=param, model=model, link=link, distr="nbinom", distrcoefs=c(size=size))$ts
  }
  fit <- tsglm(ts=timser, model=model, link=link, distr="nbinom", ...)
  result <- c(fit$coefficients, fit$distrcoefs)
  return(result)
}

distrcoef_sim <- function(N, settings){
  seeds <- sample(1e+8, size=N)
  duration <- system.time(estimates <- parSapply(cl=clust, seeds, distrcoef_sim1, size=settings$size, n=settings$n, link=settings$link))[3]
  result <- list(
    settings=settings,
    estimates=estimates,
    durations=duration
  )
  print(c(duration=duration))
  return(result)
}

N <- 1000
distrcoef_n200_size1_id <- distrcoef_sim(N=N, settings=list(n=200, link="identity", size=1))
distrcoef_n200_size1_log <- distrcoef_sim(N=N, settings=list(n=200, link="log", size=1))
distrcoef_n200_size5_id <- distrcoef_sim(N=N, settings=list(n=200, link="identity", size=5))
distrcoef_n200_size5_log <- distrcoef_sim(N=N, settings=list(n=200, link="log", size=5))
distrcoef_n200_size10_id <- distrcoef_sim(N=N, settings=list(n=200, link="identity", size=10))
distrcoef_n200_size10_log <- distrcoef_sim(N=N, settings=list(n=200, link="log", size=10))
distrcoef_n200_size20_id <- distrcoef_sim(N=N, settings=list(n=200, link="identity", size=20))
distrcoef_n200_size20_log <- distrcoef_sim(N=N, settings=list(n=200, link="log", size=20))
distrcoef_n200_sizeInf_id <- distrcoef_sim(N=N, settings=list(n=200, link="identity", size=Inf))
distrcoef_n200_sizeInf_log <- distrcoef_sim(N=N, settings=list(n=200, link="log", size=Inf))
save(
  distrcoef_n200_size1_id, distrcoef_n200_size1_log,
  distrcoef_n200_size5_id, distrcoef_n200_size5_log,
  distrcoef_n200_size10_id, distrcoef_n200_size10_log,
  distrcoef_n200_size20_id, distrcoef_n200_size20_log,
  distrcoef_n200_sizeInf_id, distrcoef_n200_sizeInf_log,
file="distrcoef_n200.RData")


N <- 200
distrcoef_n100_size1_id <- distrcoef_sim(N=N, settings=list(n=100, link="identity", size=1))
distrcoef_n100_size1_log <- distrcoef_sim(N=N, settings=list(n=100, link="log", size=1))
distrcoef_n500_size1_id <- distrcoef_sim(N=N, settings=list(n=500, link="identity", size=1))
distrcoef_n500_size1_log <- distrcoef_sim(N=N, settings=list(n=500, link="log", size=1))
distrcoef_n1000_size1_id <- distrcoef_sim(N=N, settings=list(n=1000, link="identity", size=1))
distrcoef_n1000_size1_log <- distrcoef_sim(N=N, settings=list(n=1000, link="log", size=1))
distrcoef_n2000_size1_id <- distrcoef_sim(N=N, settings=list(n=2000, link="identity", size=1))
distrcoef_n2000_size1_log <- distrcoef_sim(N=N, settings=list(n=2000, link="log", size=1))
save(
  distrcoef_n100_size1_id, distrcoef_n100_size1_log,
  distrcoef_n500_size1_id, distrcoef_n500_size1_log,
  distrcoef_n1000_size1_id, distrcoef_n1000_size1_log,
  distrcoef_n2000_size1_id, distrcoef_n2000_size1_log,
file="distrcoef_size1.RData")


# # # # #
# QIC   #
# # # # #

qic_simu <- function(n=100, link="identity"){
  timser <- tsglm.sim(n=100, param=list(intercept=ifelse(link=="identity", 4*0.5, log(4)*0.5), past_obs=0.3, past_mean=0.2), model=list(past_obs=1, past_mean=1), link=link, distr="poisson")$ts
  timser_fit <- tsglm(timser, model=list(past_obs=1, past_mean=1), link=link, distr="poisson")
  result <- c(AIC=AIC(timser_fit), QIC=QIC(timser_fit))
  return(result)
}
set.seed(0103)
system.time(aicqic_id_n100 <- data.frame(t(replicate(200, qic_simu(n=100, link="identity")))))
system.time(aicqic_log_n100 <- t(replicate(200, qic_simu(n=100, link="log"))))
save(aicqic_id_n100, aicqic_log_n100, file="qic.RData")


## # # # # # # # # # #
## Start estimation  #
## # # # # # # # # # #
#
#startest_sim1 <- function(seed, startmethod, link, n, param=list(intercept=ifelse(link=="identity", 0.5*4, 0.5*log(4)), past_obs=0.3, past_mean=0.2), ...){
#  require(tscount)
#  if(!missing(seed)) set.seed(seed)
#  model <- list(past_obs=1, past_mean=1)
#  timser <- tsglm.sim(n=n, param=param, model=model, link=link, distr="poisson")$ts
#  if(link=="identity") duration <- system.time(fit <- try(ingarch.fit(ts=timser, model=model, start.control=list(method=startmethod), ...)))[3]
#  if(link=="log") duration <- system.time(fit <- try(loglin.fit(ts=timser, model=model, start.control=list(method=startmethod), ...)))[3]
#  if(class(fit)=="try-error") return(rep(NA, 10))
#  result <- c(final=fit$coefficients, start=fit$start, loglik=fit$logLik, optimtotal=fit$final$counts[1], optimouter=fit$final$outer.iterations, duration=duration)
#  return(result)
#}
#
#startest_sim <- function(N, settings, seed=1027){
#  if(!is.null(seed)) set.seed(seed)
#  seeds <- sample(1e+8, size=N)
#  duration <- system.time(estimates <- parSapply(cl=clust, seeds, startest_sim1, startmethod=settings$startmethod, n=settings$n, link=settings$link))[3]
#  result <- list(
#    settings=settings,
#    estimates=estimates,
#    durations=duration
#  )
#  return(result)
#}
#
#N <- 200
#startest_CSS_id <- list(
#  n100=startest_sim(N=N, settings=list(startmethod="CSS", n=100, link="identity")),
#  n500=startest_sim(N=N, settings=list(startmethod="CSS", n=500, link="identity")),
#  n1000=startest_sim(N=N, settings=list(startmethod="CSS", n=1000, link="identity")),
#  n2000=startest_sim(N=N, settings=list(startmethod="CSS", n=2000, link="identity"))
#)
#startest_GLM_id <- list(
#  n100=startest_sim(N=N, settings=list(startmethod="GLM", n=100, link="identity")),
#  n500=startest_sim(N=N, settings=list(startmethod="GLM", n=500, link="identity")),
#  n1000=startest_sim(N=N, settings=list(startmethod="GLM", n=1000, link="identity")),
#  n2000=startest_sim(N=N, settings=list(startmethod="GLM", n=2000, link="identity"))
#)
#startest_iid_id <- list(
#  n100=startest_sim(N=N, settings=list(startmethod="iid", n=100, link="identity")),
#  n500=startest_sim(N=N, settings=list(startmethod="iid", n=500, link="identity")),
#  n1000=startest_sim(N=N, settings=list(startmethod="iid", n=1000, link="identity")),
#  n2000=startest_sim(N=N, settings=list(startmethod="iid", n=2000, link="identity"))
#)
#
#startest_CSS_log <- list(
#  n100=startest_sim(N=N, settings=list(startmethod="CSS", n=100, link="log")),
#  n500=startest_sim(N=N, settings=list(startmethod="CSS", n=500, link="log")),
#  n1000=startest_sim(N=N, settings=list(startmethod="CSS", n=1000, link="log")),
#  n2000=startest_sim(N=N, settings=list(startmethod="CSS", n=2000, link="log"))
#)
#startest_GLM_log <- list(
#  n100=startest_sim(N=N, settings=list(startmethod="GLM", n=100, link="log")),
#  n500=startest_sim(N=N, settings=list(startmethod="GLM", n=500, link="log")),
#  n1000=startest_sim(N=N, settings=list(startmethod="GLM", n=1000, link="log")),
#  n2000=startest_sim(N=N, settings=list(startmethod="GLM", n=2000, link="log"))
#)
#startest_iid_log <- list(
#  n100=startest_sim(N=N, settings=list(startmethod="iid", n=100, link="log")),
#  n500=startest_sim(N=N, settings=list(startmethod="iid", n=500, link="log")),
#  n1000=startest_sim(N=N, settings=list(startmethod="iid", n=1000, link="log")),
#  n2000=startest_sim(N=N, settings=list(startmethod="iid", n=2000, link="log"))
#)
#
#save(
#  startest_CSS_id, startest_CSS_log,
#  startest_GLM_id, startest_GLM_log,
#  startest_iid_id, startest_iid_log,
#file="startest.RData")
#
##load("startest.RData")
##startest_rmse_id <- list(
##  RMSE_start_CSS_id=sqrt(sapply(startest_CSS_id, function(y) apply((y$estimates[4:6,] - c(0.5*4, 0.3, 0.2))^2, 1, mean))),
##  RMSE_final_CSS_id=sqrt(sapply(startest_CSS_id, function(y) apply((y$estimates[1:3,] - c(0.5*4, 0.3, 0.2))^2, 1, mean))),
##  RMSE_start_GLM_id=sqrt(sapply(startest_GLM_id, function(y) apply((y$estimates[4:6,] - c(0.5*4, 0.3, 0.2))^2, 1, mean))),
##  RMSE_final_GLM_id=sqrt(sapply(startest_GLM_id, function(y) apply((y$estimates[1:3,] - c(0.5*4, 0.3, 0.2))^2, 1, mean))),
##  RMSE_start_iid_id=sqrt(sapply(startest_iid_id, function(y) apply((y$estimates[4:6,] - c(0.5*4, 0.3, 0.2))^2, 1, mean))),
##  RMSE_final_iid_id=sqrt(sapply(startest_iid_id, function(y) apply((y$estimates[1:3,] - c(0.5*4, 0.3, 0.2))^2, 1, mean)))
##)
##startest_rmse_log <- list(
##  RMSE_start_CSS_log=sqrt(sapply(startest_CSS_log, function(y) apply((y$estimates[4:6,] - c(0.5*log(4), 0.3, 0.2))^2, 1, mean))),
##  RMSE_final_CSS_log=sqrt(sapply(startest_CSS_log, function(y) apply((y$estimates[1:3,] - c(0.5*log(4), 0.3, 0.2))^2, 1, mean))),
##  RMSE_start_GLM_log=sqrt(sapply(startest_GLM_log, function(y) apply((y$estimates[4:6,] - c(0.5*log(4), 0.3, 0.2))^2, 1, mean))),
##  RMSE_final_GLM_log=sqrt(sapply(startest_GLM_log, function(y) apply((y$estimates[1:3,] - c(0.5*log(4), 0.3, 0.2))^2, 1, mean))),
##  RMSE_start_iid_log=sqrt(sapply(startest_iid_log, function(y) apply((y$estimates[4:6,] - c(0.5*log(4), 0.3, 0.2))^2, 1, mean))),
##  RMSE_final_iid_log=sqrt(sapply(startest_iid_log, function(y) apply((y$estimates[1:3,] - c(0.5*log(4), 0.3, 0.2))^2, 1, mean)))
##)
##
##startest_rmseplot <- function(x, index, label, ylim=c(0.3, 1.3)){
##  samplesizes <- c(100, 500, 1000, 2000)
##  plot(NA, xlim=c(100, 2000), ylim=ylim, type="n", xlab="Length of time series", ylab="RMSE", log="x", xaxt="n", cex.lab=1.4)
##  lines(samplesizes, x[[1]][index,], type="o", lty="dashed")
##  lines(samplesizes, x[[2]][index,], type="o")
##  lines(samplesizes, x[[3]][index,], type="o", pch=2, lty="dashed")
##  lines(samplesizes, x[[4]][index,], type="o", pch=2)
##  lines(samplesizes, x[[5]][index,], type="o", pch=8, lty="dashed")
##  lines(samplesizes, x[[6]][index,], type="o", pch=8)
##  text(120, 0, label=label, cex=2, pos=3)
##  #legend("bottomleft", legend="", title=label, bty="n", cex=2)
##}
##
##par(mfrow=c(3,1), mar=c(0.25,4,0,0), las=1, mgp=c(2.5,0.6,0), oma=c(2.5,0,2.5,1))
##startest_rmseplot(startest_rmse_id, index=1, label=expression(beta[0]), ylim=c(0, 1))
##startest_rmseplot(startest_rmse_id, index=2, label=expression(beta[1]), ylim=c(0, 0.3))
##startest_rmseplot(startest_rmse_id, index=3, label=expression(alpha[1]), ylim=c(0, 0.25))
##axis(side=1, line=0)
##mtext(text="Length of time series", side=1, line=1.7, cex=0.9)
##title(main="         Linear model", outer=TRUE, cex.main=1.6)
##par(mfrow=c(3,1), mar=c(0.25,4,0,0), las=1, mgp=c(2.5,0.6,0), oma=c(2.5,0,2.5,1))
##startest_rmseplot(startest_rmse_log, index=1, label=expression(beta[0]), ylim=c(0, 1))
##startest_rmseplot(startest_rmse_log, index=2, label=expression(beta[1]), ylim=c(0, 0.3))
##legend("topright", legend=c("ARMA ", "GLM ", "ARMA", "GLM"), ncol=2, title="Initial estim.:   Final estim.:    ", lwd=1, lty=c("dashed", "dashed", "solid", "solid"), pch=c(1,2,1,2), bg="white", cex=1.1)
##startest_rmseplot(startest_rmse_log, index=3, label=expression(alpha[1]), ylim=c(0, 0.4))
##axis(side=1, line=0)
##mtext(text="Length of time series", side=1, line=1.7, cex=0.9)
##title(main="        Log-linear model", outer=TRUE, cex.main=1.6)
###Caption: Root mean square error (RMSE) of two different start estimators (dashed lines) for a linear (left) respectively log-linear (right) model of order $p=q=1$. The time series are simulated from the same model. The final maximum likelihood estimators based on the respective start estimator are shown with solid lines. Each value is based on 1000 replications.



## # # # # # # # # # # # # # # # #
## Comparison with GLM and ARMA  #
## # # # # # # # # # # # # # # # #
#
##Is it worth using tsglm or can the functions glm (ignoring the feedback) or arima (ignoring the integer-valued nature of the data) obtain good results as well?
#
#comparison_sim1 <- function(seed, link, n, param=list(intercept=ifelse(link=="identity", 0.1*1, 0.1*log(1)), past_obs=0.5, past_mean=0.4), ...){
#  shortpoisint <- function(lambda, conf_level=0.95){
#    dens <- dpois(0:100, lambda)
#    set <- (order(dens)[cumsum(sort(dens)) >= 1-conf_level]) - 1
#    interval <- range(set)
#    return(interval)
#  }
#  require(tscount)
#  if(!missing(seed)) set.seed(seed)
#  model <- list(past_obs=1, past_mean=1)
#  timser <- tsglm.sim(n=n+1, param=param, model=model, link=link, distr="poisson")$ts
#  true <- timser[n+1]
#  timser <- timser[1:n]
#  fit_our <- try(tsglm(ts=timser, model=model, link=link, distr="poisson", ...))
#  if(class(fit_our)=="try-error"){
#    pred_our <- NA
#    predint_our <- c(NA, NA)
#  }else{
#    pred_our <- predict(fit_our)$pred[1]
#    predint_our <- shortpoisint(pred_our) #highest density region interval
#    #predint_our <- pred_our + c(-1,+1)*qnorm(0.975)*sqrt(pred_our) #normal approximation
#  }
#  dat <- data.frame(obs=timser[-1], past_obs=timser[-n])
#  fit_GLM <- try(glm(obs ~ 1 + past_obs, data=dat, family=poisson(link=link), start=c(1,0)))
#  if(any(class(fit_GLM)=="try-error")){
#    pred_GLM <- NA
#    predint_GLM <- c(NA, NA)
#  }else{
#    pred_GLM <- predict(fit_GLM, newdata=data.frame(past_obs=timser[n]), se.fit=TRUE, type="response")$fit[[1]]
#    predint_GLM <- shortpoisint(pred_GLM) #highest density region interval
#    #predint_GLM <- pred_GLM + c(-1,+1)*qnorm(0.975)*sqrt(pred_GLM) #normal approximation
#  }
#  fit_ARMA <- try(arima(timser, order=c(1, 0, 1), method="ML"))
#  if(class(fit_ARMA)=="try-error"){
#    pred_ARMA <- NA
#    predint_ARMA <- c(NA, NA)
#  }else{
#    pred_ARMA <- predict(fit_ARMA)
#    predint_ARMA <- pred_ARMA$pred[1] + c(-1, +1)*qnorm(0.975)*pred_ARMA$se[1]
#    predint_ARMA[1] <- max(predint_ARMA[1], 0)
#    predint_ARMA[2] <- ceiling(predint_ARMA[2])
#  }
#  result <- c(true=true, our=pred_our, our=predint_our, GLM=pred_GLM, GLM=predint_GLM, ARMA=pred_ARMA$pred[1], ARMA=predint_ARMA)
#  return(result)
#}
#
#comparison_sim <- function(N, settings){
#  seeds <- sample(1e+8, size=N)
#  duration <- system.time(values <- t(parSapply(cl=clust, seeds, comparison_sim1, n=settings$n, link=settings$link)))[3]
#  result <- list(
#    settings=settings,
#    values=values,
#    durations=duration
#  )
#  print(c(duration=duration))
#  return(result)
#}
#
#N <- 1000
#comparison_n200_id <- comparison_sim(N=N, settings=list(n=200, link="identity"))
#comparison_n200_log <- comparison_sim(N=N, settings=list(n=200, link="log"))
#
#save(
#  comparison_n200_id,
#  comparison_n200_log,
#file="comparison.RData")
#
##load("comparison.RData")
##comp <- comparison_n200_id$values
###comp <- comparison_n200_log$values
####Root mean square prediction error:
##sqrt(mean((comp[, "true"] - comp[, "our"])^2, na.rm=TRUE))
##sqrt(mean((comp[, "true"] - comp[, "GLM"])^2, na.rm=TRUE))
##sqrt(mean((comp[, "true"] - comp[, "ARMA"])^2, na.rm=TRUE))
####Coverage rate of the prediction interval:
##mean(comp[, "true"] >= comp[, "our1"] & comp[, "true"] <= comp[, "our2"], na.rm=TRUE)
##mean(comp[, "true"] >= comp[, "GLM1"] & comp[, "true"] <= comp[, "GLM2"], na.rm=TRUE)
##mean(comp[, "true"] >= comp[, "ARMA1"] & comp[, "true"] <= comp[, "ARMA2"], na.rm=TRUE)
####Length of the prediction interval:
##mean(comp[, "our2"] - comp[, "our1"], na.rm=TRUE)
##mean(comp[, "GLM2"] - comp[, "GLM1"], na.rm=TRUE)
##mean(comp[, "ARMA2"] - comp[, "ARMA1"], na.rm=TRUE)
###For both, the linear and the log-linear model our findings are as follows: glm has a poor RMSPE as compared to tsglm, hence it does not predict the mean very well because it lacks the feedback. arima has to wide prediction intervals because it does not use the appropriate distribution.


######################################################################

sessionInfo()
Sys.time()
memory.size()
gc()

#Shut down cluster for parallel computation:
stopCluster(clust)
mpi.exit()

