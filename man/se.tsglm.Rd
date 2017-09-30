\name{se.tsglm}
\alias{se.tsglm}
\alias{se}
\title{
Standard Errors of a Fitted Generalised Linear Model for Time Series of Counts
}
\description{
Computes the standard errors for the parameters of a fitted GLM-type model for time series of counts.
}
\usage{
\method{se}{tsglm}(object, B, parallel, level=0.95, ...)
}
\arguments{
  \item{object}{
an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}.
}
  \item{B}{
positive integer value giving the number of bootstrap samples to use for estimation of the standard errors. If missing the standard errors are based on a normal approximation.
  }
  \item{parallel}{
  logical value. If \code{parallel=TRUE}, the bootstrap is distributed to multiple cores parallely. Requires a computing cluster to be initialised and registered as the default cluster by \code{\link[parallel:makeCluster]{makeCluster}} and \code{\link[parallel:setDefaultCluster]{setDefaultCluster}} from package \code{parallel}. 
}
  \item{level}{
numeric value determining the desired coverage rate of confidence intervals.
}
  \item{...}{
  additional arguments to be passed to the fitting function \code{\link{tsglm}}. Only made use of if the standard errors are computed by a bootstrap procedure.
}
}

\details{
By default the standard errors and confidence intervals are based on a normal approximation of the (quasi) maximum likelihood estimator. The standard errors are the square roots of the diagonal elements of the inverse of the information matrix. Because there is no analytical approximation of the standard error for the overdispersion coefficient \code{sigmasq}, its standard error and its confidence interval are set to \code{NA}.

If the number of bootstrap samples \code{B} is given, the standard errors and condidence intervals are computed by a parametric bootstrap. The standard errors are the empirical standard deviation of the parameter estimations of \code{B} random samples drawn from the fitted model given in argument \code{object}. The confidence intervals are the \code{a}- and \code{(1-a)}-quantile of this bootstrap sample with \code{a=(1-level)/2}.
}

\value{
A list with the following components:

\item{est}{a vector of the maximum likelihood estimated coefficients.}

\item{se}{a vector of the standard errors of each estimated coefficient.}

\item{ci}{a matrix with the columns \code{"lower"} and \code{"upper"} giving the lower and upper boundaries of confidence intervals for the model parameters.}

\item{level}{numerical value giving the coverage rate of the confidence intervals.}

\item{type}{a character value \code{"normapprox"} or \code{"bootstrap"} giving how the standard errors are computed.}

If the standard errors are computed by a parametric bootstrap procedure, the following component is additionally returned:

\item{B}{positive integer value giving the number of bootstrap samples used for estimation of the standard errors.}
}

\references{
Liboschik, T., Fokianos, K. and Fried, R. (201?) tscount: An R package for analysis of count time series following generalized linear models. \emph{Journal of Statistical Software}, to appear. Available as a vignette: \code{vignette("tsglm", package="tscount")}.
}

\author{
Tobias Liboschik and Philipp Probst
}

\seealso{
\code{\link{tsglm}} for fitting a GLM for time series of counts.
}

\examples{
###Road casualties in Great Britain (see help("Seatbelts"))
timeseries <- Seatbelts[, "VanKilled"]
regressors <- cbind(PetrolPrice=Seatbelts[, c("PetrolPrice")],
                    linearTrend=seq(along=timeseries)/12)
#Logarithmic link function with Poisson distribution:                    
seatbeltsfit <- tsglm(ts=timeseries, link="log",
  model=list(past_obs=c(1, 12)), xreg=regressors, distr="poisson")
  
se(seatbeltsfit) #by normal approximation

\dontrun{
system.time(stderror <- se(seatbeltsfit, B=100)) #by bootstrap
stderror
#This estimation of bootstrap standard errors takes several minutes on a single
#processing unit, of course depending on its speed.
#Parallel computation for shorter run time on a cluster:
library(parallel)
ntasks <- 3
clust <- makeCluster(ntasks)
setDefaultCluster(cl=clust)
system.time(stderror <- se(seatbeltsfit, B=100, parallel=TRUE))}
}

\keyword{Inference}
