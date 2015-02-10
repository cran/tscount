\name{ingarch.analytical}
\alias{ingarch.analytical}
\alias{ingarch.mean}
\alias{ingarch.acf}
\alias{ingarch.var}
\title{
Analytical Mean, Variance and Autocorrelation of an INGARCH Process 
}
\description{
Functions to calculate the analytical mean, variance and autocorrelation / partial autocorrelation / autocovariance function of an integer-valued generalised autoregressive conditional heteroscedasticity (INGARCH) process.
}
\usage{
ingarch.mean(intercept, past_obs=NULL, past_mean=NULL)
ingarch.var(intercept, past_obs=NULL, past_mean=NULL)
ingarch.acf(intercept, past_obs=NULL, past_mean=NULL, lag.max=10,
        type=c("acf", "pacf", "acvf"), plot=TRUE, ...)
}

\arguments{
  \item{intercept}{
  numeric positive value for the intercept \eqn{\beta_0}{\beta[0]}.
}
  \item{past_obs}{
  numeric non-negative vector containing the coefficients \eqn{\beta_1,\ldots, \beta_p}{\beta[1], \ldots, \beta[p]} for regression on previous observations (see Details).
}
  \item{past_mean}{
  numeric non-negative vector containing the coefficients \eqn{\alpha_1,\ldots, \alpha_q}{\alpha[1], \ldots, \alpha[q]} for regression on previous conditional means (see Details).
}
  \item{lag.max}{
  integer value indicating how many lags of the (partial) autocorrelation / autocovariance function should be calculated.
}
  \item{type}{
  character. If \code{type="acf"} (the default) the autocorrelation function is calculated, \code{"pacf"} gives the partial autocorrelation function and \code{"acvf"} the autocovariance function.
}
  \item{plot}{
  logical. If \code{plot=TRUE} (the default) the values are plotted and returned invisible. 
}
  \item{...}{
  additional arguments to be passed to function \code{\link{plot}}.
  }
}

\details{
  The INGARCH model of order \eqn{p} and \eqn{q} used here follows the definition
    \deqn{Z_{t}|{\cal{F}}_{t-1} \sim \mathrm{Poi}(\kappa_{t}),}{
    Z[t]|F[t-1] ~ Poi(\kappa[t]),}
  where \eqn{{\cal{F}}_{t-1}}{F[t-1]} is the history of the process up to time \eqn{t-1} and \eqn{\mathrm{Poi}}{Poi} is the Poisson distribution parametrised by its mean (cf. Ferland et al., 2006).
  The conditional mean \eqn{\kappa_t}{\kappa[t]} is given by
    \deqn{\kappa_t = \beta_0 + \beta_1 Z_{t-1} + \ldots + \beta_p Z_{t-p}
    + \alpha_1 \kappa_{t-1} + \ldots  + \alpha_q \kappa_{t-q}.}{
    \kappa[t] = \beta[0] + \beta[1] Z[t-1] + \ldots + \beta[p] Z[t-p]
    + \alpha[1] \kappa[t-1] + \ldots  + \alpha[q] \kappa[t-q].}
    The function \code{ingarch.acf} depends on the function \code{\link[ltsa:tacvfARMA]{tacvfARMA}} from package \code{ltsa}, which needs to be installed.
}

\author{
Tobias Liboschik
}

\references{
Ferland, R., Latour, A. and Oraichi, D. (2006) Integer-valued GARCH process. \emph{Journal of Time Series Analysis} \bold{27(6)}, 923--942, \url{http://dx.doi.org/10.1111/j.1467-9892.2006.00496.x}.
}

\seealso{
\code{\link{tsglm}} for fitting a more genereal GLM for time series of counts of which this INGARCH model is a special case. \code{\link{tsglm.sim}} for simulation from such a model.
}

\examples{
ingarch.mean(0.3, c(0.1,0.1), 0.1)
\dontrun{
ingarch.var(0.3, c(0.1,0.1), 0.1)
ingarch.acf(0.3, c(0.1,0.1,0.1), 0.1, type="acf", lag.max=15)}
}
