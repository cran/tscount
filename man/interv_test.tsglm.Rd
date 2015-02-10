\name{interv_test.tsglm}
\alias{interv_test.tsglm}
\alias{interv_test}
\title{
Testing for Interventions in Count Time Series Following Generalised Linear Models
}

\description{
Test for one or more interventions of given type at given time as proposed by Fokianos and Fried (2010, 2012).
}

\usage{
\method{interv_test}{tsglm}(fit, tau, delta, external, 
            info=c("score"), est_interv=FALSE, ...)
}

\arguments{
  \item{fit}{
an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}.
}
  \item{tau}{
  integer vector of times at which the interventions occur which are tested for.
}
  \item{delta}{
  numeric vector that determines the types of the interventions (see Details). Must be of the same length as \code{tau}.
}
  \item{external}{
  logical vector of length \code{length(tau)} specifying for each intervention wether its effect is external or not (see Details). If this is only a scalar this choice will be used for all interventions. If this is only a scalar this choice will be used for all interventions. If omitted all interventions will have an internal effect (i.e. \code{external=FALSE}).
} 
  \item{info}{
  character value that determines how to calculate the information matrix, see \code{\link{tsglm}}. Currently \code{"score"} is the only possible choice.
}
  \item{est_interv}{
  logical value. If \code{est_interv=TRUE} a fit for the model with all specified interventions is computed and additionally returned.
}
  \item{...}{
  additional arguments passed to the fitting function \code{\link{tsglm}}.
}
}
 
\details{
A score test on the null hypothesis of no interventions is done. The null hypothesis is that the data are generated from the model specified in the argument \code{model}, see definition in \code{\link{tsglm}}. Under the alternative there are one or more intervention effects occuring at times \code{tau}. The types of the intervention effects are specified by \code{delta} as defined in \code{\link{interv_covariate}}. The interventions are included as additional covariates according to the definition in \code{\link{tsglm}}. It can have an internal (the default) or external (\code{external=TRUE}) effect (see Liboschik et al., 2014).

Under the null hypothesis the test statistic has asymptotically a chi-square distribution with \code{length(tau)} (i.e. the number of breaks) degrees of freedom. The returned p-value is based on this and approximately valid for long time series, i.e. when \code{length(ts)} large.
}

\value{
An object of class \code{"interv_test"}, which is a list with at least the following components:

  \item{test_statistic}{
  value of the test statistic.
}
  \item{df}{
  degrees of freedom of the chi-squared distribution the test statistic is compared with.
}
  \item{p-value}{
  p-value of the test.
}
  \item{fit_H0}{
  object of class \code{"tsglm"} with the fitted model under the null hypothesis of no intervention, see \code{\link{tsglm}}.
}
  \item{model_interv}{
  model specification of the model with the specified interventions.
}

If argument \code{est_interv=TRUE}, the following component is additionally returned:

  \item{fit_interv}{
  object of class \code{"tsglm"} with the fitted model with the specified interventions, see \code{\link{tsglm}}.
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
S3 method \code{\link{print}}.

\code{\link{tsglm}} for fitting a GLM for time series of counts.
\code{\link[=interv_detect.tsglm]{interv_detect}} for detection of single interventions of given type and \code{\link[=interv_multiple.tsglm]{interv_multiple}} for iterative detection of multiple interventions of unknown types. \code{\link{interv_covariate}} for generation of deterministic covariates describing intervention effects.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
#Test for the intervention effects which were found in Fokianos und Fried (2010):
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
campyfit_intervtest <- interv_test(fit=campyfit, tau=c(84,100), delta=c(1,0))
campyfit_intervtest
}

\keyword{Intervention detection}
