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
\method{se}{tsglm}(object, B, parallel, ...)
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
  \item{...}{
  additional arguments to be passed to the fitting function \code{\link{tsglm}}. Only made use of if the standard errors are computed by a bootstrap procedure.
}
}

\details{
By default the standard errors are based on a normal approximation of the (quasi) maximum likelihood estimator. They are the square roots of the diagonal elements of the inverse of the information matrix. There is no analytical approximation of the standard error for the overdispersion coefficient \code{sigmasq}, which is therefore set to \code{NA}.

If the number of bootstrap samples \code{B} is given, the standard errors are computed by a parametric bootstrap. The standard errors are the empirical standard deviation of the parameter estimations of \code{B} random samples drawn from the fitted model given in argument \code{object}.
}

\value{
A list with the following components:

\item{est}{a vector of the maximum likelihood estimated coefficients.}

\item{se}{a vector of the standard errors of each estimated coefficient.}

\item{type}{a character value \code{"normapprox"} or \code{"bootstrap"} giving how the standard errors are computed.}

If the standard errors are computed by a parametric bootstrap procedure, the following component is additionally returned:

\item{B}{positive integer value giving the number of bootstrap samples used for estimation of the standard errors.}
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
  model=list(past_obs=c(1, 12)), xreg=regressors, distr="pois")

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
