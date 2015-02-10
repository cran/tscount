\name{predict.tsglm}
\alias{predict.tsglm}
\title{
Predicts Method for Time Series of Counts Following Generalised Linear Models
}
\description{
Predict future observations based on a fitted GLM-type model for time series of counts.
}
\usage{
\method{predict}{tsglm}(object, n.ahead = 1, newobs = NULL, newxreg = NULL, level=0.95, B, ...)
}
\arguments{
  \item{object}{
an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}.
}
  \item{n.ahead}{
positive integer value giving the number of steps ahead for which predictions should be made.
}
  \item{newobs}{
integer vector of known future observations of the time series. This argument is only relevant if more than one observation ahead is to be predicted (\code{n.ahead} greater than 1). The \eqn{h}-step-ahead prediction for \eqn{h>1} is computed as a 1-step-ahead prediction given all previous values, which can be observations of the original time series or new observations provided in this argument. Previous observations which are not available are replaced by their respective 1-step-ahead prediction.
}
  \item{newxreg}{
matrix or vector containing new values for the covariates to be used for prediction. If \code{newxreg} is omitted or contains less rows than the value of \code{n.ahead}, the last known values of the covariates are used for prediction. This is often not reasonable and it is strongly advised to explicitely make assumptions on future covariates and to specify the argument \code{xreg} accordingly. 
}
  \item{level}{
numeric value determining the confidence level for prediction intervals. Note that this level holds only for each of the prediction intervals but not globally.
}
  \item{B}{
positive integer value giving the number of bootstrap samples to use for numerical determination of prediction intervals with a confidence level of \code{level} each. If missing no prediction intervals are computed.
}
  \item{...}{
further arguments are currently ignored. Only for compatibility with generic function.
}
}
\details{
Returns predictions for the \code{n.ahead} observations following the fitted time series contained in argument \code{object}. The 1-step-ahead prediction is the conditional expectation of the observation to be predicted given the past. The true parameters are replaced by their estimations given in argument \code{object}. For a 2-step-ahead-prediction the true previous observation is used when given in argument \code{newobs}, otherwise it is replaced by the 1-step-ahead prediction computed before. For a 3-step-prediction this holds for the previous two observations, which are replaced by their respective predictions if not available, and so on.

The function additionally returns prediction intervals numerically determined by a parametric bootstrap, if the number of bootstrap replicates is specified by the argument \code{B}. Note that the prediction intervals do not reflect the additional uncertainty induced by the parameter estimation. However, for sufficiently long time series it is expected that this uncertainty is negligible compared to the uncertainty of the predictive distribution.
}
\value{
A list with the following element:

  \item{pred}{a numeric vector of the predictions. Has class \code{"ts"} if the response used for fitting has this class.}

If the predictive distribution is approximated by simulation (i.e. argument \code{B} specifying the number of simulated trajectories of the process is provided), the list has the additional elements:

  \item{median}{a vector giving the empirical median of the simulated predictive distribution for each of the future time points. Has class \code{"ts"} if the response used for fitting has this class.}

  \item{interval_shortest}{a matrix with the columns \code{"lower"} and \code{"upper"} giving the lower and upper boundaries of prediction intervals for the future time points, each with an intended coverage rate as given in argument \code{level}. The prediction intervals are obtained from the simulated predictive distribution such that they have minimal length. Has class \code{"ts"} if the response used for fitting has this class.}
  
  \item{interval_quantiles}{like \code{interval_shortest}, but the prediction intervals are given by the empirical \code{a}- and (1-\code{a})-quantiles of the simulated predictive distribution, with \code{a}=(1-\code{level})/2.}

  \item{type}{a character value saying which method was used for computing the prediction intervals. Currently this is always the method \code{"bootstrap"}.}

  \item{B}{an integer value giving the number of bootstrap samples which were used for computing prediction intervals.}
}

\author{
Tobias Liboschik and Philipp Probst
}

\seealso{
\code{\link{tsglm}} for fitting a GLM for time series of counts.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
predict(campyfit, n.ahead=5)
}

\keyword{Prediction}
