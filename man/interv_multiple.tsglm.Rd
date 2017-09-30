\name{interv_multiple.tsglm}
\alias{interv_multiple.tsglm}
\alias{interv_multiple}

\title{
Detecting Multiple Interventions in Count Time Series Following Generalised Linear Models
}
\description{
Iterative detection procedure for multiple interventions of unknown types occuring at unknown times as proposed by Fokianos and Fried (2010, 2012).
}
\usage{
\method{interv_multiple}{tsglm}(fit, taus=2:length(fit$ts), deltas=c(0,0.8,1),
                external=FALSE, B=10, signif_level=0.05,
                start.control_bootstrap, final.control_bootstrap,
                inter.control_bootstrap, parallel=FALSE, ...)
}

\arguments{
  \item{fit}{
an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}.
}
  \item{taus}{
  integer vector of times which are considered for the possible intervention to occur.
  Default is to consider all times.
}
  \item{deltas}{
  numeric vector that determines the types of intervention to be considered (see Details).
}
  \item{external}{
  logical value specifying wether the interventions effect is external or not (see Details).
} 
  \item{B}{
  positive integer value giving the number of bootstrap samples for estimation of the p-value.
}
  \item{signif_level}{
  numeric value with \code{0 <= signif_level <= 1} giving a significance level for the procedure.
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
  \item{...}{
  additional arguments passed to the function for detection of single intervention effects \code{\link[=interv_detect.tsglm]{interv_detect}} and via this function some of the arguments are passed to the fitting function \code{\link{tsglm}}.
}
}

\details{
This function performs an iterative procedure for detection of multiple intervention effects. In each step the function \code{\link[=interv_detect.tsglm]{interv_detect}} is applied for each of the possible intervention types provided in the argument \code{deltas}. If there is (after a Bonferroni correction) no significant intervention effect the procedure stops. Otherwise the type of intervention with the minimum p-value is chosen. In case of equal p-values preference is given to a level shift (i.e. \eqn{\delta=1}) and then to the type of intervention with the largest test statistic. The effect of the chosen intervention is removed from the time series. The time series cleaned from the intervention effect is tested for further interventions in a next step.

For each time in \code{taus} the test statistic of a score test on an intervention effect occuring at that time is computed, see \code{\link[=interv_test.tsglm]{interv_test}}. The time with the maximum test statistic is considered as a candidate for a possible intervention effect at that time. The type of the intervention effect is specified by \code{delta} as described in \code{\link{interv_covariate}}. The intervention is included as an additional covariate according to the definition in \code{\link{tsglm}}. It can have an internal (the default) or external (\code{external=TRUE}) effect (see Liboschik et al., 2014).

All p-values given in the output are multiplied by the number of intervention types considered to account for the multiple testing in each step by a Bonferroni correction. Note that this correction can lead to p-values greater than one.

Note that this bootstrap procedure is very time-consuming.
}

\value{
An object of class \code{"interv_multiple"}, which is a list with the following components:

  \item{interventions}{
  data frame giving the detected interventions, which has the variables \code{tau}, \code{delta}, \code{size}, \code{test_statistic} and \code{p-value}. 
}
  \item{fit_H0}{
  object of class \code{"tsglm"} with the fitted model under the null hypothesis of no intervention, see \code{\link{tsglm}}.
}
  \item{fit_cleaned}{
  object of class \code{"tsglm"} with the fitted model for the cleanded time series after the last step of the iterative procedure, see \code{\link{tsglm}}.
}
  \item{model_interv}{
  model specification of the model with all detected interventions at their respective times.
}
  \item{fit_interv}{
  object of class \code{"tsglm"} with the fitted model with all detected interventions at their respective times, see \code{\link{tsglm}}.
}
  \item{track}{
  named list of matrices with the detailed results of the iterative detection procedure. Element \code{tau_max} gives the times where the test statistic has its maximum for each type of intervention and in each iteration step and element \code{size} gives the estimated sizes of the respective intervention effects. Elements \code{test_statistic} and \code{p_value} require no further explanation.
}
}

\references{
Fokianos, K. and Fried, R. (2010) Interventions in INGARCH processes. \emph{Journal of Time Series Analysis} \bold{31(3)}, 210--225, \url{http://dx.doi.org/10.1111/j.1467-9892.2010.00657.x}.

Fokianos, K., and Fried, R. (2012) Interventions in log-linear Poisson autoregression. \emph{Statistical Modelling} \bold{12(4)}, 299--322. \url{http://dx.doi.org/10.1177/1471082X1201200401}.

Liboschik, T. (2016) Modelling count time series following generalized linear models. \emph{PhD Thesis TU Dortmund University}, \url{http://dx.doi.org/10.17877/DE290R-17191}.

Liboschik, T., Kerschke, P., Fokianos, K. and Fried, R. (2016) Modelling interventions in INGARCH processes. \emph{International Journal of Computer Mathematics} \bold{93(4)}, 640--657, \url{http://dx.doi.org/10.1080/00207160.2014.949250}.
}

\author{
Tobias Liboschik, Philipp Probst, Konstantinos Fokianos and Roland Fried
}

\seealso{
S3 methods \code{\link{print}} and \code{\link[=plot.interv_multiple]{plot}}.

\code{\link{tsglm}} for fitting a GLM for time series of counts.
\code{\link[=interv_test.tsglm]{interv_test}} for testing for intervention effects and \code{\link[=interv_detect.tsglm]{interv_detect}} for detection of single interventions of given type. \code{\link{interv_covariate}} for generation of deterministic covariates describing intervention effects.
}

\examples{
\dontrun{
###Campylobacter infections in Canada (see help("campy"))
#Searching for potential intervention effects (runs several hours!):
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
campyfit_intervmultiple <- interv_multiple(fit=campyfit, taus=80:120,
                              deltas=c(0,0.8,1), B=500, signif_level=0.05)
campyfit_intervmultiple
plot(campyfir_intervmultiple)
#Parallel computation for shorter run time on a cluster:
library(parallel)
ntasks <- 3
clust <- makeCluster(ntasks)
setDefaultCluster(cl=clust)
interv_multiple(fit=campyfit, taus=80:120, deltas=c(0,0.8,1), B=500,
                signif_level=0.05, parallel=TRUE)}
}

\keyword{Intervention detection}
