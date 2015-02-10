\name{interv_detect.tsglm}
\alias{interv_detect.tsglm}
\alias{interv_detect}

\title{
Detecting an Intervention in Count Time Series Following Generalised Linear Models
}

\description{
Detection procedure for an intervention of given type occuring at unknown time as proposed by Fokianos and Fried (2010, 2012).
}

\usage{
\method{interv_detect}{tsglm}(fit, taus=2:length(fit$ts), delta, external=FALSE, 
              B=NULL, info=c("score"), start.control_bootstrap, 
              final.control_bootstrap, inter.control_bootstrap, 
              parallel=FALSE, est_interv=TRUE, ...)
}

\arguments{
  \item{fit}{
an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}.
}
  \item{taus}{
  integer vector of time points which are considered for the possible intervention to occur.
  Default is to consider all possible time points.
}
  \item{delta}{
  numeric value that determines the type of intervention (see Details).
}
  \item{external}{
  logical value specifying wether the intervention's effect is external or not (see Details).
} 
  \item{B}{
  positive integer value giving the number of bootstrap samples for estimation of the p-value. For \code{B=NULL} (the default) no p-value is returned.
}
  \item{info}{
  character value that determines how to calculate the information matrix, see \code{\link{tsglm}}. Currently \code{"score"} is the only possible choice.
}
  \item{start.control_bootstrap}{
  named list that determines how to make initial estimation in the bootstrap, see argument \code{start.control} in \code{\link{tsglm}}. If missing, the same settings as for the regular estimation are used.
}
  \item{final.control_bootstrap}{
  named list that determines how to make final maximum likelihood estimation in the bootstrap, see argument \code{final.control} in \code{\link{tsglm}}. If missing, the same settings as for the regular estimation are used. If \code{final.control_bootstrap=NULL}, then the model is not re-fitted for each bootstrap sample. Instead the parameters of the original fit which have been used for simulating the bootstrap samples are used. This approach saves computation time at the cost of a more conservative procedure, see Fokianos and Fried (2012).
}
  \item{inter.control_bootstrap}{
  named list  determining how to maximise the log-likelihood function in an intermediate step, see argument \code{inter.control} in \code{\link{tsglm}}. If missing, the same settings as for the regular estimation are used.
}
  \item{parallel}{
  logical value. If \code{parallel=TRUE}, the bootstrap is distributed to multiple cores parallely. Requires a computing cluster to be initialised and registered as the default cluster by \code{\link[parallel:makeCluster]{makeCluster}} and \code{\link[parallel:setDefaultCluster]{setDefaultCluster}} from package \code{parallel}.
}
  \item{est_interv}{
  logical value. If \code{est_interv=TRUE} a fit for the model with the intervention effect with the largest test statistic is computed and additionally returned.
}
  \item{...}{
  additional arguments passed to the fitting function \code{\link{tsglm}}.
}
}
 
\details{
For each time in \code{taus} the score test statistic for an intervention effect occuring at that time is computed, see \code{\link[=interv_test.tsglm]{interv_test}}. The time with the maximum test statistic is considered as a candidate for a possible intervention effect at that time. The type of the intervention effect is specified by \code{delta} as described in \code{\link{interv_covariate}}.  The intervention is included as an additional covariate according to the definition in \code{\link{tsglm}}. It can have an internal (the default) or external (\code{external=TRUE}) effect (see Liboschik et al., 2014).

If argument \code{B} is not \code{NULL}, the null hypothesis that there is no intervention effect at any time is tested. Test statistic for this test is the maximum test statistic of the score test (see above). The p-value is computed by a parametric bootstrap with \code{B} bootstrap samples. It is recommended to use at least several hundred bootstrap samples. Note that this bootstrap procedure is very time-consuming.
}

\value{
  An object of class \code{"interv_detect"}, which is a list with at least the following components:

  \item{test_statistic}{
  maximum value of the score test statistics for all considered times in \code{taus}.
}
  \item{test_statistic_tau}{
  numeric vector of all score test statistics at the considered times in \code{taus}.
}
  \item{tau_max}{
  time at which the score test statistic has its maximum.
}
  \item{fit_H0}{
  object of class \code{"tsglm"} with the fitted model under the null hypothesis of no intervention, see \code{\link{tsglm}}.
}
  \item{model_interv}{
  model specification of the model with the specified intervention at time \code{tau_max}.
}

If argument \code{est_interv=TRUE} (the default), the following component is additionally returned:

  \item{fit_interv}{
  object of class \code{"tsglm"} with the fitted model with the specified intervention at time \code{tau_max}, see \code{\link{tsglm}}.
}
}

\references{
Fokianos, K. and Fried, R. (2010) Interventions in INGARCH processes. \emph{Journal of Time Series Analysis} \bold{31(3)}, 210--225, \url{http://dx.doi.org/10.1111/j.1467-9892.2010.00657.x}.

Fokianos, K., and Fried, R. (2012) Interventions in log-linear Poisson autoregression. \emph{Statistical Modelling} \bold{12(4)}, 299--322. \url{http://dx.doi.org/10.1177/1471082X1201200401}.

Liboschik, T., Kerschke, P., Fokianos, K. and Fried, R. (2014) Modelling interventions in INGARCH processes. \emph{International Journal of Computer Mathematics} (published online), \url{http://dx.doi.org/10.1080/00207160.2014.949250}.
}

\author{
Tobias Liboschik, Philipp Probst, Konstantinos Fokianos and Roland Fried
}

\seealso{
S3 methods \code{\link{print}} and \code{\link[=plot.interv_detect]{plot}}.

\code{\link{tsglm}} for fitting a GLM for time series of counts.
\code{\link[=interv_test.tsglm]{interv_test}} for testing on intervention effects and \code{\link[=interv_multiple.tsglm]{interv_multiple}} for iterative detection of multiple interventions of unknown types. \code{\link{interv_covariate}} for generation of deterministic covariates describing intervention effects.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
#Searching for a potential intervention effect:
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
campyfit_intervdetect <- interv_detect(fit=campyfit, taus=80:120, delta=1)
campyfit_intervdetect
plot(campyfit_intervdetect)
#Additionally computing a p-value with the bootstrap procedure based on 500
#replications would take about 20 minutes in this example on a single
#processing unit, of course depending on its speed.

\dontrun{
#Parallel computation for shorter run time on a cluster:
library(parallel)
ntasks <- 3
clust <- makeCluster(ntasks)
setDefaultCluster(cl=clust)
interv_detect(fit=campyfit, taus=80:120, delta=1, B=500, parallel=TRUE)}
}

\keyword{Intervention detection}
