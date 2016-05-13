\name{pit}
\alias{pit}
\alias{pit.default}
\alias{pit.tsglm}
\title{
Predictive Model Assessment with a Probability Integral Transform Histogram
}
\description{
The function allows a probabilistic calibration check with a Probability Integral Transform (PIT) histogram.
}
\usage{
\method{pit}{tsglm}(object, bins=10, ...)
\method{pit}{default}(response, pred, distr=c("poisson", "nbinom"), distrcoefs, bins=10, ...)
}
\arguments{
  \item{object}{
an object of class \code{"tsglm"}.
}
  \item{bins}{
number of bins in the histogram. Default value is 10.
}
  \item{response}{
integer vector. Vector of observed values.
}
  \item{pred}{
numeric vector. Vector of predicted values.
}
  \item{distr}{
character giving the conditional distribution. Currently implemented are the Poisson (\code{"poisson"})and the Negative Binomial (\code{"nbinom"}) distribution.
}
  \item{distrcoefs}{
numeric vector of additional coefficients specifying the conditional distribution. For \code{distr="poisson"} no additional parameters need to be provided. For \code{distr="nbinom"} the additional parameter \code{size} needs to be specified (e.g. by \code{distrcoefs=2}), see \code{\link{tsglm}} for details.
}
%  \item{ci}{
%coverage probability for confidence band. Plotting of the confidence interval is suppressed if \code{ci} is not between zero and one.
%}
%  \item{ci.col}{
%colour to plot the confidence band lines.
%}
  \item{...}{
additional arguments passed to \code{\link{plot}}.
}
}
\details{
A PIT histogram is a tool for evaluating the statistical consistency between the probabilistic forecast and the observation. The predictive distributions of the observations are compared with the actual observations. If the predictive distribution is ideal the result should be a flat PIT histogram with no bin having an extraordinary high or low level. For more information about PIT histograms see the references listed below.

%The condidence band plotted is based on the assumption of an ideal forecaster, under which the PIT values are independently and identically uniformly distributed on \eqn{[0,1]}{[0,1]}. It holds the given coverage probability for all bars simultaneously.
}

\references{
Christou, V. and Fokianos, K. (2013) On count time series prediction. \emph{Journal of Statistical Computation and Simulation} (published online), \url{http://dx.doi.org/10.1080/00949655.2013.823612}.

Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment for count data. \emph{Biometrics} \bold{65}, 1254--1261, \url{http://dx.doi.org/10.1111/j.1541-0420.2009.01191.x}.

Gneiting, T., Balabdaoui, F. and Raftery, A.E. (2007) Probabilistic forecasts, calibration and sharpness. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} \bold{69}, 243--268, \url{http://dx.doi.org/10.1111/j.1467-9868.2007.00587.x}.
}

\author{
Philipp Probst and Tobias Liboschik
}

\seealso{
\code{\link{tsglm}} for fitting a GLM for time series of counts.

\code{\link{marcal}} and \code{\link{scoring}} for other predictive model assessment tools.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
pit(campyfit)
}

\keyword{Model assessment}
