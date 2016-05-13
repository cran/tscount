\name{scoring}
\alias{scoring}
\alias{scoring.default}
\alias{scoring.tsglm}
\title{
Predictive Model Assessment with Proper Scoring Rules
}
\description{
Computes scores for the assessment of sharpness of a fitted model for time series of counts.
}
\usage{
\method{scoring}{tsglm}(object, individual=FALSE, cutoff=1000, ...)
\method{scoring}{default}(response, pred, distr=c("poisson", "nbinom"), distrcoefs,
          individual=FALSE, cutoff=1000, ...)
}
\arguments{
  \item{object}{
an object of class \code{"tsglm"}.
}
  \item{individual}{
logical. If \code{FALSE} (the default) the average scores are returned. Otherwise a matrix with the individual scores for each observation is returned.
}
  \item{cutoff}{
positive integer. Summation over the infinite sample space \{0,1,2,...\} of a distribution is cut off at this value. This affects the quadratic, spherical and ranked probability score.
}
  \item{response}{
integer vector. Vector of observed values \eqn{Y_1,...,Y_n}{Y[1],...,Y[n]}.
}
  \item{pred}{
numeric vector. Vector of predicted values \eqn{\mu_{P_1},...,\mu_{P_n}}{\mu_P[1],...,\mu_P[n]}.
}
  \item{distr}{
character giving the conditional distribution. Currently implemented are the Poisson (\code{"poisson"})and the Negative Binomial (\code{"nbinom"}) distribution.
}
  \item{distrcoefs}{
numeric vector of additional coefficients specifying the conditional distribution. For \code{distr="poisson"} no additional parameters need to be provided. For \code{distr="nbinom"} the additional parameter \code{size} needs to be specified (e.g. by \code{distrcoefs=2}), see \code{\link{tsglm}} for details.
}
  \item{...}{
further arguments are currently ignored. Only for compatibility with generic function.
}
}

\details{
The scoring rules are penalties that should be minimised for a better forecast, so a smaller scoring value means better sharpness. 
Different competing forecast models can be ranked via these scoring rules. 
They are computed as follows:
For each score \eqn{s} and time \eqn{t} the value \eqn{s(P_{t},Y_{t})}{s(P[t],Y[t])} is computed, where \eqn{P_t}{P[t]} is the predictive 
c.d.f. and \eqn{Y_t}{Y[t]} is the observation at time \eqn{t}. To obtain the overall score for one model the average of the score of all observations
\eqn{(1/n) \sum_{t=1}^{n}s(P_{t},Y_{t})}{(1/n) \sum s(P[t],Y[t])}
is calculated. 

For all \eqn{t \geq 1}, let \eqn{p_{y} = P(Y_{t}=y | {\cal{F}}_{t-1} )}{p[y]=P(Y[t]=y | F[t-1])} be the density function of the predictive distribution at \eqn{y} and 
\eqn{||p||^2=\sum_{y=0}^{\infty} p_y^2}{||p||^2= \sum p[y]^2} be a quadratic sum over the whole sample space \eqn{y=0,1,2,...} of the predictive distribution.
\eqn{\mu_{P_t}}{\mu_P[t]} and \eqn{\sigma_{P_t}}{\sigma_P[t]} are the mean and the standard deviation of the predictive distribution, respectively.

Then the scores are defined as follows:

Logarithmic score: \eqn{logs(P_{t},Y_{t})= -log p_{y}}{logs(P[t],Y[t])= -log p[y] }

Quadratic or Brier score: \eqn{qs(P_{t},Y_{t}) = -2p_{y} + ||p||^2 }{qs(P[t],Y[t])= -2p[y] + ||p||^2}

Spherical score: \eqn{sphs(P_{t},Y_{t})=\frac{-p_{y}}{||p||}}{sphs(P[t],Y[t])= -p[y] / ||p||}

Ranked probability score:  \eqn{rps(P_{t},Y_{t})=\sum_{x=0}^{\infty}(P_{t}(x) - 1(Y_t\leq x))^2}{rps(P[t],Y[t])=\sum (P[t](x) - 1(Y[t]\le x))^2} (sum over the whole sample space \eqn{x=0,1,2,...})

Dawid-Sebastiani score: \eqn{dss(P_{t},Y_{t})=\left(\frac{Y_t-\mu_{P_t}}{\sigma_{P_t}}\right)^2 + 2log\sigma_{P_t}}{dss(P[t],Y[t]) = ( (Y[t]-\mu_P[t]) / (\sigma_P[t]) )^2 + 2 log \sigma_P[t]}

Normalized squared error score: \eqn{nses(P_{t},Y_{t})=\left(\frac{Y_t-\mu_{P_t}}{\sigma_{P_t}}\right)^2}{nses(P[t],Y[t])= ( (Y[t]-\mu_P[t]) \ (\sigma_P[t]) )^2}

Squared error score: \eqn{ses(P_{t},Y_{t})=(Y_t-\mu_{P_t})^2}{ses(P[t],Y[t])=(Y[t]-\mu_P[t])^2}

For more information on scoring rules see the references listed below.

}

\value{ 
Returns a named vector of the mean scores (if argument \code{individual=FALSE}, the default) or a data frame of the individual scores for each observation (if argument \code{individual=TRUE}). The scoring rules are named as follows:

\item{logarithmic}{Logarithmic score}

\item{quadratic}{Quadratic or Brier score}

\item{spherical}{Spherical score}

\item{rankprob}{Ranked probability score}

\item{dawseb}{Dawid-Sebastiani score}

\item{normsq}{Normalized squared error score}

\item{sqerror}{Squared error score}
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

\code{\link{pit}} and \code{\link{marcal}} for other predictive model assessment tools.

\code{\link[surveillance]{permutationTest}} in package \code{surveillance} for the Monte Carlo permutation test for paired individual scores by Paul and Held (2011, \emph{Statistics in Medicine} \bold{30}, 1118--1136, \url{http://dx.doi.org/10.1002/sim.4177}).
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
scoring(campyfit)
}

\keyword{Model assessment}
