\name{tsglm.sim}
\alias{tsglm.sim}
\title{
Simulate a Time Series Following a Generalised Linear Model
}
\description{
Generates a simulated time series from a GLM-type model for time series of counts (see \code{\link{tsglm}} for details).
}
\usage{
tsglm.sim(n, param = list(intercept = 1, past_obs = NULL, past_mean = NULL,
            xreg = NULL), model = list(past_obs = NULL, past_mean = NULL,
            external = FALSE), xreg = NULL, link = c("identity", "log"),
            distr = c("poisson", "nbinom"), distrcoefs, fit, n_start = 50)
}
\arguments{
  \item{n}{
  integer value giving the number of observations to be simulated.
}
  \item{param}{
  a named list giving the parameters for the linear predictor of the model, which has the following elements: 
    \describe{
      \item{\code{intercept}}{numeric positive value for the intercept \eqn{\beta_0}{\beta[0]}.}
      \item{\code{past_obs}}{numeric non-negative vector containing the coefficients \eqn{\beta_1,\ldots, \beta_p}{\beta[1], \ldots, \beta[p]} for regression on previous observations (see Details).}
      \item{\code{past_mean}}{numeric non-negative vector containing the coefficients \eqn{\alpha_1,\ldots, \alpha_q}{\alpha[1], \ldots, \alpha[q]} for regression on previous conditional means (see Details).}
      \item{\code{xreg}}{numeric non-negative vector specifying the size \eqn{\nu_1,\ldots, \nu_r}{\nu[1], \ldots, \nu[r]} of each intervention}   
  }
}
  \item{model}{
  a named list specifying the model for the linear predictor, which has the elements \code{past_obs}, \code{past_mean} and \code{external} (see function \code{\link{tsglm}} for details). This model specification must be in accordance to the parameters given in argument \code{param}.
}
  \item{xreg}{
  matrix with covariates in the columns (see \code{\link{tsglm}} for details). Its number of rows must be equal to the number of observations which should be simulated.
}
  \item{link}{
    character giving the link function. Default is \code{"identity"}, simulating from a so-called INGARCH model. Another possible choice is \code{"log"}, simulating from a log-linear model.
}
  \item{distr}{
    character giving the conditional distribution. Default is \code{"poisson"}, i.e. a Poisson distribution.
}
  \item{distrcoefs}{
  numeric vector of additional coefficients specifying the conditional distribution. For \code{distr="poisson"} no additional parameters need to be provided. For \code{distr="nbinom"} the additional parameter \code{size} needs to be specified (i.e. by \code{distrcoefs=2}), see \code{\link{tsglm}} for details.
}
  \item{fit}{
  an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}. If argument \code{fit} is not missing, the specification of the linear predictor, the link function and the estimated parameters from this argument are used instead of those in arguments \code{model}, \code{link} and \code{param}. The length of the simulated time series is only taken from argument \code{fit}, if no argument \code{n} is provided. The same holds for arguments \code{xreg}, \code{distr} and \code{distrcoefs}, which are also prefered over the respective information provided in argument \code{fit} if both are provided.
}
  \item{n_start}{
  number of observations used as a burn-in.
}
}

\details{
The definition of the model used here is like in function \code{\link{tsglm}}.

Note that during the burn-in period covariates are set to zero.

If a previous model fit is given in argument \code{fit} and the length of the burn-in period \code{n_start} is set to zero, then the a continuation of the original time series is simulated.
}

\value{
  A list with the following components:

\item{\code{ts}}{
  an object of class \code{"ts"} with the simulated time series.
}
\item{\code{linear.predictors}}{
  an object of class \code{"ts"} with the simulated conditional means \eqn{\kappa_t}{\kappa[t]} for all \eqn{t=1, \ldots, n}.
}
\item{\code{xreg.effects}}{
  an object of class \code{"ts"} with the cumulated effect of the covariates \eqn{\eta_1 X_{t,1} + \ldots + \eta_r X_{t,r}}{\eta[1] X[t,1] + \ldots + \eta[r] X[t,r]} for all \eqn{t=1, \ldots, n}.
}
}

\author{
Tobias Liboschik and Philipp Probst
}

\seealso{
\code{\link{tsglm}} for fitting a GLM for time series of counts.
}

\examples{
#Simulate from an INGARCH model with two interventions:
interventions <- interv_covariate(n=200, tau=c(50, 150), delta=c(1, 0.8))
model <- list(past_obs=1, past_mean=c(1, 7), external=FALSE)
param <- list(intercept=2, past_obs=0.3, past_mean=c(0.2, 0.1), xreg=c(3, 10))
tsglm.sim(n=200, param=param, model=model, xreg=interventions, link="identity",
          distr="nbinom", distrcoefs=c(size=1))
}

\keyword{Simulation}
