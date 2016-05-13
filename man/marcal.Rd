\name{marcal}
\alias{marcal}
\alias{marcal.default}
\alias{marcal.tsglm}
\title{
Predictive Model Assessment with a Marginal Calibration Plot
}
\description{
The function produces a marginal calibration plot.
}
\usage{
\method{marcal}{tsglm}(object, plot=TRUE, ...)
\method{marcal}{default}(response, pred, distr=c("poisson", "nbinom"), distrcoefs, plot=TRUE, ...)
}
\arguments{
  \item{object}{
an object of class \code{"tsglm"}.
}
  \item{plot}{
  logical. If \code{plot=TRUE} (the default), the marginal calibration is plotted and the underlying data are returned invisibly only. 
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
  \item{...}{
additional arguments to be passed to \code{\link{plot}}.
}
}
\details{
Marginal Calibration can be assessed by taking the difference between the average predictive cumulative distribution function (c.d.f.) and the empirical c.d.f. of the observations. Minor fluctuations about zero are expected if the marginal calibration hypothesis is true. For more information about marginal calibration see the refererences listed below.
}
\value{
Produces a plot of the difference between the average predictive cumulative distribution function (c.d.f.) and the empirical c.d.f. of the observations at each value between the highest and lowest observation of the time series (only for \code{plot=TRUE}).

Returns a list with elements \code{x} and \code{y}, where \code{x} are the threshold values and \code{y} the respective differences of predictive and empirical cumulative distribution function (invisibly for \code{plot=TRUE}).
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

\code{\link{pit}} and \code{\link{scoring}} for other predictive model assessment tools.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
marcal(campyfit)
}

\keyword{Model assessment}
