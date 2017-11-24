\name{predict.tsglm}
\alias{predict.tsglm}
\title{
Predicts Method for Time Series of Counts Following Generalised Linear Models
}
\description{
Predict future observations based on a fitted GLM-type model for time series of counts.
}
\usage{
\method{predict}{tsglm}(object, n.ahead=1, newobs=NULL, newxreg=NULL, level=0.95,
        global=FALSE, type=c("quantiles", "shortest", "onesided"),
        method=c("conddistr", "bootstrap"), B=1000,
        estim=c("ignore", "bootstrap", "normapprox", "given"), B_estim=B,
        coefs_given, ...)
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
matrix or vector containing new values for the covariates to be used for prediction. If \code{newxreg} is omitted or contains less rows than the value of \code{n.ahead}, the last known values of the covariates are used for prediction. This is usually not reasonable and it is strongly advised to explicitely make assumptions on future covariates and to specify the argument \code{xreg} accordingly. 
}
  \item{level}{
numeric value determining the desired coverage rate of prediction intervals. If \code{level=0} no prediction intervals are computed.
}
  \item{global}{
logical value saying whether the coverage rate for \eqn{Y_{n+1},...,Y_{n+h}}{Y[n+1],...,Y[n+h]} specified by argument \code{level} holds globally (\code{global=TRUE}) or for each of the \code{n.ahead} prediction intervals individually (\code{global=FALSE}, the default). In the former case the individual coverage rate for a single prediction interval is Bonferroni adjusted to a level of \code{1-(1-level)/n.ahead}.
}
  \item{type}{
character value saying how the prediction interval shall be constructed. If \code{type="quantiles"} (the default), its limits are chosen to be the \code{a}- and \code{(1-a)}-quantiles of the respective (approximated) distribution, with \code{a=(1-level)/2}. If \code{type="shortest"} it is chosen such that it has minimal length. Note that these two types of construction principles frequently lead to the same result. If \code{type="onesided"} a one-sided prediction interval is constructed where the lower boundary is always zero.
}
  \item{method}{
character value saying which method to be used for computing the prediction intervals. If \code{method="conddistr"} the prediction intervals are based on the conditional distribution given by the model with the unknown parameters being replaced by their respective estimations. This is only possible if only 1-step-ahead predictions are to be computed (possibly recursively using the new observations given in argument \code{newobs}). If \code{method="bootstrap"} the predictive distribution is approximated by a parametric bootstrap where \code{B} trajectories of the process are simulated from the fitted model. This is currently only possible if no new observations are given in argument \code{newobs}. By default the method \code{"conddistr"} is prefered whenever it is appliccable.
}
  \item{B}{
positive integer value giving the number of samples of a parametric bootstrap to use for numerical determination of prediction intervals (only necessary if argument \code{method="bootstrap"}).
}
  \item{estim}{
character value saying how the prediction intervals shall account for the additional uncertainty induced by the parameter estimation. This is particularly important if the model was fitted on a short time series. If \code{estim="ignore"} (the default), this additional uncertainty is ignored. The other two options (\code{estim="bootstrap"} and \code{estim="normapprox"}) are only possible if \code{method="bootstrap"}. If these are selected the bootstrap samples are not generated from a model with the parameters of the original fit. Instead, each of the \code{B} bootstrap samples is generated from a model with parameters which are itself randomly generated. This two-stage approach should take into account the additional estimation uncertainty.
If \code{estim="bootstrap"}, the parameters are obtained from a fit to a parametric bootstrap replication of the original time series.
If \code{estim="normapprox"}, the regression parameters are generated from a multivariate normal distribution which is based on the normal approximation of the original quasi maximum likelihood estimator and reflects the estimation uncertainty. In that case the additional distribution coefficients are not randomly generated such that their estimation uncertainty is ignored.
If \code{estim="given"}, the parameters are resampled from a table of possible parameters which need to be given in argument \code{coefs_given}.
}
  \item{B_estim}{
positive integer value giving the number of parameters used for resampling to account for estimation uncertainty. Only necessary for \code{estim="bootstrap"} and \code{estim="normapprox"}. If \code{B_estim} is smaller than \code{B}, the parameters are resampled with replacement.
}
  \item{coefs_given}{
table with parameters in the rows. Only necessary for \code{estim="given"}). If \code{nrow(coefs_given)} is smaller than \code{B}, the parameters are resampled with replacement.
}
  \item{...}{
further arguments are currently ignored. Only for compatibility with generic function.
}
}
\details{
Returns predictions for the \code{n.ahead} observations following the fitted time series contained in argument \code{object}. The 1-step-ahead prediction is the conditional expectation of the observation to be predicted given the past. The true parameters are replaced by their estimations given in argument \code{object}. For a 2-step-ahead-prediction the true previous observation is used when given in argument \code{newobs}, otherwise it is replaced by the 1-step-ahead prediction computed before. For a 3-step-prediction this holds for the previous two observations, which are replaced by their respective predictions if not available, and so on.

Unless \code{level=0}, the function also returns prediction intervals. Read the description of the arguments \code{type} and\code{method} for further details on the computation. Note that the prediction intervals do not reflect the additional uncertainty induced by the parameter estimation. However, for sufficiently long time series used for model fitting, it is expected that this uncertainty is negligible compared to the uncertainty of the predictive distribution. The argument \code{estim} allows to account fot this additional estimation uncertainty if \code{method="bootstrap"}, see the description of this argument.

If prediction intervals are computed the function additionally returns the median of the predictive distribution. If \code{method="conddistr"} this is the analytical median of the conditional distribution, otherwise the empirical median of the simulated distribution.
}
\value{
A list with at least the following element:

  \item{pred}{a numeric vector of the predictions. Has class \code{"ts"} if the response used for fitting has this class.}
  
If prediction intervals are calculated, the list has the additional element:

  \item{interval}{a matrix with the columns \code{"lower"} and \code{"upper"} giving the lower and upper boundaries of prediction intervals for the future time points, each with an intended coverage rate as given in argument \code{level}. Has class \code{"ts"} if the response used for fitting has this class.}
  
  \item{level}{a numeric value determining the desired coverage rate of prediction intervals.}
  
  \item{global}{a logical value saying whether the coverage rate \code{level} holds globally or for each of the prediction intervals individually.}

  \item{type}{a character value saying how the prediction intervals were computed. Possible values are \code{"quantiles"} and \code{"shortest"}.}
  
  \item{method}{a character value saying which method were used for computation of prediction intervals. Possible values are \code{"conddistr"} and \code{"bootstrap"}.}
  
  \item{B}{an integer value giving the number of bootstrap samples which were used for computing prediction intervals. Is \code{NULL} if computation was done by \code{method="conddistr"}.}
  
  \item{estim}{a character value saying how the prediction intervals account for estimation uncertainty of the model parameters. Possible values are \code{"ignore"}, \code{"bootstrap"}, \code{"normapprox"} and \code{"given"}.}
  
  \item{B_estim}{an integer value giving the number of parameter values used for resampling to account for estimation uncertainty. This value is zero if the estimation uncertainty is ignored.}
  
  \item{warning_messages}{a character vector containing warning messages. This should be \code{NULL} if no warning messages occured.}
  
  \item{median}{a vector giving the median of the predictive distribution for each of the future time points. Has class \code{"ts"} if the response used for fitting has this class.}
}

\references{
Liboschik, T., Fokianos, K. and Fried, R. (2017) tscount: An R package for analysis of count time series following generalized linear models. \emph{Journal of Statistical Software} \bold{82(5)}, 1--51, \url{http://dx.doi.org/10.18637/jss.v082.i05}.
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
predict(campyfit, n.ahead=1) #prediction interval using conditional distribution
predict(campyfit, n.ahead=5, global=TRUE) #prediction intervals using parametric bootstrap
}

\keyword{Prediction}
