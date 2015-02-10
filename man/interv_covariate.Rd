\name{interv_covariate}
\alias{interv_covariate}

\title{
Describing Intervention Effects for Time Series with Deterministic Covariates
}
\description{
Generates covariates describing certain types of intervention effects according to the definition by Fokianos and Fried (2010).
}
\usage{
interv_covariate(n, tau, delta)
}
\arguments{
  \item{n}{
integer value giving the number of observations the covariates should have.
}
  \item{tau}{
integer vector giving the times where intervention effects occur.
}
  \item{delta}{
numeric vector with constants specifying the type of intervention (see Details). Must be of the same length as \code{tau}.
}
}

\details{
The intervention effect occuring at time \eqn{\tau} is described by the covariate
\deqn{X_t = \delta^{t-\tau} I_{[\tau,\infty)}(t),}{X_t = \delta^(t-\tau) I(t>=\tau),}
where \eqn{I_{[\tau,\infty)}(t)}{I(t>=\tau)} is the indicator function which is 0 for \eqn{t < \tau} and 1 for \eqn{t \geq \tau}{t >= \tau}. The constant \eqn{\delta} with \eqn{0 \leq \delta \leq 1}{0 <= \delta <= 1} specifies the type of intervention. For \eqn{\delta = 0} the intervention has an effect only at the time of its occurence, for \eqn{0 < \delta < 1} the effect decays exponentially and for \eqn{\delta = 1} there is a persistent effect of the intervention after its occurence.

If \code{tau} and \code{delta} are vectors, one covariate is generated with \code{tau[1]} as \eqn{\tau} and \code{delta[1]} as \eqn{\delta}, another covariate for the second elements and so on.
}

\value{
A matrix with \code{n} rows and \code{length(tau)} columns. The generated covariates describing the interventions are the columns of the matrix.
}

\references{
Fokianos, K. and Fried, R. (2010) Interventions in INGARCH processes. \emph{Journal of Time Series Analysis} \bold{31(3)}, 210--225, \url{http://dx.doi.org/10.1111/j.1467-9892.2010.00657.x}.

Fokianos, K., and Fried, R. (2012) Interventions in log-linear Poisson autoregression. \emph{Statistical Modelling} \bold{12(4)}, 299--322. \url{http://dx.doi.org/10.1177/1471082X1201200401}.

Liboschik, T., Kerschke, P., Fokianos, K. and Fried, R. (2014) Modelling interventions in INGARCH processes. \emph{International Journal of Computer Mathematics} (published online), \url{http://dx.doi.org/10.1080/00207160.2014.949250}.
}

\author{
Tobias Liboschik
}

\seealso{
\code{\link{tsglm}} for fitting a GLM for time series of counts.
\code{\link[=interv_test.tsglm]{interv_test}}, \code{\link[=interv_detect.tsglm]{interv_detect}} and \code{\link[=interv_multiple.tsglm]{interv_multiple}} for tests and detection procedures for intervention effects.
}

\examples{
interv_covariate(n=140, tau=c(84,100), delta=c(1,0))
}

\keyword{Intervention detection}
