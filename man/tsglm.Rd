\name{tsglm}
\alias{tsglm}
\alias{print.tsglm}
\alias{logLik.tsglm}
\alias{tsglm.meanfit}
\alias{vcov.tsglm}
\title{
Count Time Series Following Generalised Linear Models
}
\description{
The function \code{tsglm} fits a generalised linear model (GLM) for time series of counts.
The specification of the linear predictor allows for regressing on past observations, past values of the linear predictor and covariates as defined in the Details section.
There is the so-called INGARCH model with the identity link (see for example Ferland et al., 2006, Fokianos et al., 2009) and another model with the logarithmic link (see for example Fokianos and Tjostheim, 2011), which also differ in the specification of the linear predictor.
The conditional distribution can be chosen to be either Poisson or negative binomial.

Estimation is done by conditional maximum likelihood for the Poisson distribution or by a conditional quasi-likelihood approach based on the Poisson likelihood function for the negative binomial distribution.

There is a vignette available which introduces the functionality of \code{tsglm} and related functions of this package and its underlying statistical methods (\code{vignette("tsglm", package="tscount")}).

The function \code{tsglm.meanfit} is a lower level function to fit the mean specification of such a model assuming a Poisson distribution. It is called by \code{tsglm}. It has additional arguments allowing for a finer control of the fitting procedure, which can be handed over from the function \code{tsglm} by its \code{...} argument. Note that it is usually not necessary for a user to call this lower level functions nor to worry about the additional arguments provided by this function. The defaults of these arguments have been chosen wisely by the authors of this package and should perform well in most applications.
}
\usage{
tsglm(ts, model = list(past_obs = NULL, past_mean = NULL,
            external = NULL), xreg = NULL, link = c("identity", "log"),
            distr = c("poisson", "nbinom"), ...)

tsglm.meanfit(ts, model, xreg, link, score = TRUE,
            info = c("score", "none", "hessian", "sandwich"),
            init.method=c("marginal", "iid", "firstobs", "zero"),
            init.drop = FALSE, epsilon = 1e-06, slackvar = 1e-06,
            start.control = list(), final.control = list(),
            inter.control = NULL)
}


\arguments{
  \item{ts}{
    a univariate time series.
  }
  \item{model}{
    a named list specifying the model for the linear predictor, which can be of the following elements:
    \describe{  
      \item{\code{past_obs}}{integer vector giving the previous observations to be regressed on (autoregression). This is a vector with the elements \eqn{i_1,\ldots,i_p}{i[1],\ldots,i[p]} (see Details). If omitted, or of length zero, there will be no regression on previous observations.}    
      \item{\code{past_mean}}{integer vector giving the previous conditional means to be regressed on. This is a vector with the elements \eqn{j_1,\ldots,j_q}{j[1],\ldots,j[q]} (see Details). If omitted, or of length zero, there will be no regression on previous conditional means.} 
      \item{\code{external}}{logical vector of length \code{ncol(xreg)} specifying for each covariate wether its effect should be external or not (see Details). If this is a scalar this choice will be used for all covariates. If omitted, all covariates will have an internal effect (i.e. \code{external=FALSE}).}
    } 
  }
  \item{xreg}{
    matrix with covariates in the columns, i.e. its number of rows must be \code{length(ts)}. This is the matrix \eqn{X} (see Details). If omitted no covariates will be included. For the identity link the covariates have to be non-negative.
  }
  \item{link}{
    character giving the link function. Default is \code{"identity"}, fitting an INGARCH model. Another possible choice is \code{"log"}, fitting a log-linear model.
  }
  \item{distr}{
    character giving the conditional distribution. Default is \code{"poisson"}, i.e. a Poisson distribution.
  }
  \item{...}{
    additional arguments to be passed to the lower level fitting function \code{tsglm.meanfit}. See below.  
  }
  \item{score}{
    logical value indicating whether the score vector should be computed.
  }
  \item{info}{
    character that determines if and how to compute the information matrix. Can be set to \code{"score"} (the default) for calculation via the outer product of the score vector, or to \code{"hessian"} for calculation via the Hessian matrix of second derivatives. For \code{info="sandwich"} the information matrix is estimated by a sandwich formula using both the outer score product and the Hessian matrix. If set to \code{"none"}, no information matrix is computed. For \code{distr="nbinom"} one can only use \code{info="score"}.
  }
  \item{init.method}{
    character that determines how the recursion of the conditional mean (and possibly of its derivatives) is initialised. If set to \code{"marginal"} (the default), the marginal mean of a model without covariates and its derivatives are used. If set to \code{"iid"}, all values are initialised by the marginal mean under the assumption of i.i.d. data, which depends on the intercept only. If set to \code{"firstobs"} the first obersvation is used. If set to \code{"zero"}, the recursions are initialised by the value zero.
  }
  \item{init.drop}{
    logical value that determines which observations are considered for computation of the log-likelihood, the score vector and, if applicable, the information matrix. If \code{TRUE}, the first \code{max(model$past_obs)} observations, which are needed for the autoregression, are not considered. If \code{FALSE} (the default), all observations are considered and pre-sample values determined by the method specified by the argument \code{itit.method} are used for the autoregression. Note that in the first case the effective number of observations used for maximum likelihood estimation is lower than the total number of observations of the original time series. Consequently only this lower number of observations is considered in the output. Note that for \code{init.drop=TRUE} the log-likelihood function for models of different orders might not be comparable if the effective number of observations is different.
  }
  \item{epsilon}{
    numeric positive but small value determining how close the parameters may come to the limits of the parameter space.
  }
  \item{slackvar}{
    numeric positive but small value determining how true inequalities among the parameter restrictions are treated; a true inequality \code{x < y} will be transformed to \code{x + slackvar <= y}.
  }
  \item{start.control}{
    named list with optional elements that determine how to make the start estimation. Possible list elements are:
    \describe{  
      \item{\code{use}}{integer vector of length one or two giving the number of observations from the beginning (if of length one) or the range of observations (if of length two) used for start estimation. For \code{use = Inf} all observations are used, which is the default.}
      \item{\code{method}}{character specifying how start estimators should be estimated. Possible values are \code{"iid"}, \code{"CSS"}, \code{"CSS-ML"}, \code{"ML"}, \code{"MM"}, \code{"GLM"} and \code{"fixed"}. 
  If \code{method} is \code{"iid"} (the default), a moment estimator assuming an iid model without covariates is used.
  If \code{method="MM"}, the start estimate is an ARMA(1,1) fit by moment estimators and parameters of higher order than one are set to zero. For this method the starting parameter values for the covariates are zero by default and can be set by the list element \code{xreg}.
  If \code{method} is \code{"CSS"}, \code{"CSS-ML"} or \code{"ML"}, the start estimate is based on an ARMA fit using the function \code{\link{arima}}, and list element \code{method} is passed to its argument of the same name.  
  If \code{method="GLM"}, the estimated parameters of a generalised linear model with regression on the specified past observations and covariates, but not on past conditional means, are used as start estimates. Initial estimates for the coefficients of past conditional means are set to zero.
  If \code{method="fixed"}, parameters given in further named list elements of \code{start.control} are used when available, else the predefined values given in the following are used.}
      \item{\code{intercept}}{numeric value with the start value for the intercept parameter. Default value is 1.}
      \item{\code{past_obs}}{numeric vector with the start values for parameters for regression on previous observations. Default values are zero.}
      \item{\code{past_mean}}{numeric vector with the start values for parameters for regression on previous conditional means. Default values are zero.}
      \item{\code{xreg}}{numeric vector with the start values for the regression parameters. These values will also be used if \code{method="MM"}. Default values are zero.}
    }  
  }
  \item{final.control}{
    named list with optional elements that determine how to make the final maximum likelihood estimation. If \code{final.control=NULL}, only start estimates are computed and a list with fewer elements which has not the class \code{"tsglm"} is returned. Possible list elements of this argument are:
    \describe{  
      \item{\code{constrained}}{named list whose elements are passed to function \code{\link{constrOptim}} with possible elements \code{mu}, \code{outer.iterations} and \code{outer.eps} (see \code{\link{constrOptim}} for details). If \code{constrained=NULL}, an unconstrained optimisation is made with function \code{\link{optim}}. Note that this is likely to result in a fitted model which is non-stationary, which might cause further problems.}
      \item{\code{optim.method}}{character which is passed to functions \code{\link{constrOptim}} or \code{\link{optim}} as argument \code{method}. The default is \code{"BFGS"}.}
      \item{\code{optim.control}}{named list which is passed to function \code{\link{constrOptim}} or \code{\link{optim}} as the argument \code{control}. Must not contain the list element \code{fnscale}. The default is \code{list(maxit=20, reltol=1e-8)}.}
    }
  }
  \item{inter.control}{
    named list determining how to maximise the log-likelihood function in a first step. This intermediate optimisation will start from the start estimation and be followed by the final optimisation, which will in turn start from the intermediate optimisation result. This intermediate optimisation is intended to use a very quick but imprecise optimisation algorithm. Possible elements are the same as for \code{final.control}. The default is \code{inter.control=NULL}, which skips this intermediate optimisation step. 
  }
}

\details{
The INGARCH model (argument \code{link="identity"}) used here follows the definition
    \deqn{Z_{t}|{\cal{F}}_{t-1} \sim \mathrm{Poi}(\nu_{t}) \quad \mathrm{or} \quad Z_{t}|{\cal{F}}_{t-1} \sim \mathrm{NegBin}(\nu_{t}, \phi),}{
    Z[t]|F[t-1] ~ Poi(\nu[t]) or  Z[t]|F[t-1] ~ NegBin(\nu[t], \phi),}
where \eqn{{\cal{F}}_{t-1}}{F[t-1]} denotes the history of the process up to time \eqn{t-1}, \eqn{\mathrm{Poi}}{Poi} and \eqn{\mathrm{NegBin}}{NegBin} is the Poisson respectively the negative binomial distribution with the parametrisation as specified below.
For the model with covariates having an internal effect (the default) the linear predictor of the INGARCH model (which is in that case identical to the conditional mean) is given by
    \deqn{\nu_t = \beta_0 + \beta_1 Z_{t-i_1} + \ldots + \beta_p Z_{t-i_p}
    + \alpha_1 \nu_{t-j_1} + \ldots  + \alpha_q \nu_{t-j_q}
    + \eta_1 X_{t,1} + \ldots + \eta_r X_{t,r}.}{
    \nu[t] = \beta[0] + \beta[1] Z[t-i[1]] + \ldots + \beta[p] Z[t-i[p]]
    + \alpha[1] \nu[t-j[1]] + \ldots  + \alpha[q] \nu[t-j[q]]
    + \eta[1] X[t,1] + \ldots + \eta[r] X[t,r].}

The log-linear model (argument \code{link="log"}) used here follows the definition
    \deqn{Z_{t}|{\cal{F}}_{t-1} \sim \mathrm{Poi}(\lambda_{t}) \quad \mathrm{or} \quad Z_{t}|{\cal{F}}_{t-1} \sim \mathrm{NegBin}(\lambda_{t}, \phi),}{
    Z[t]|F[t-1] ~ Poi(\lambda[t]) or  Z[t]|F[t-1] ~ NegBin(\lambda[t], \phi),}
with \eqn{\lambda_t = \exp(\nu_t)}{\lambda[t] = \exp(\nu[t])} and \eqn{{\cal{F}}_{t-1}}{F[t-1]} as above.
For the model with covariates having an internal effect (the default) the linear predictor \eqn{\nu_t = \log(\lambda_t)}{\nu[t] = \log(\lambda[t])} of the log-linear model is given by
    \deqn{\nu_t = \beta_0 + \beta_1 \log(Z_{t-i_1}+1) + \ldots + \beta_p \log(Z_{t-i_p}+1)
    + \alpha_1 \nu_{t-j_1} + \ldots  + \alpha_q \nu_{t-j_q}
    + \eta_1 X_{t,1} + \ldots + \eta_r X_{t,r}.}{
    \nu[t] = \beta[0] + \beta[1] \log(Z[t-i[1]]+1) + \ldots + \beta[p] \log(Z[t-i[p]]+1)
    + \alpha[1] \nu[t-j[1]] + \ldots  + \alpha[q] \nu[t-j[q]]
    + \eta[1] X[t,1] + \ldots + \eta[r] X[t,r].}
Note that because of the logarithmic link function the effect of single summands in the linear predictor on the conditional mean is multiplicative and hence the parameters play a different role than in the INGARCH model, although they are denoted by the same letters.
  
The Poisson distribution is parametrised by the mean \code{lambda} according to the definition in \code{\link{Poisson}}.
The negative binomial distribution is parametrised by the mean \code{mu} with an additional dispersion parameter \code{size} according to the definition in \code{\link{NegBinomial}}. In the notation above its mean parameter \code{mu} is \eqn{\nu_{t}}{\nu[t]} and its dispersion parameter \code{size} is \eqn{\phi}.

This function allows to include covariates in two different ways. A covariate can have a so-called internal effect as defined above, where its effect propagates via the regression on past values of the linear predictor and on past observations. Alternatively, it can have a so-called external effect, where its effect does not directly propagates via the feedback on past values of the linear predictor, but only via past observations.
For external effects of the covariates, the linear predictor for the model with identity link is given by
    \deqn{\nu_t = \mu_t + \eta_1 X_{t,1} + \ldots + \eta_r X_{t,r},}{
    \nu[t] = \mu[t] + \eta[1] X[t,1] + \ldots + \eta[r] X[t,r],}
    \deqn{\mu_t = \beta_0 + \beta_1 Z_{t-i_1} + \ldots + \beta_p Z_{t-i_p}
    + \alpha_1 \mu{t-j_1} + \ldots  + \alpha_q \mu{t-j_q},}{
    \mu[t] = \beta[0] + \beta[1] Z[t-i[1]] + \ldots + \beta[p] Z[t-i[p]
    + \alpha[1] \mu[t-j[1]] + \ldots  + \alpha[q] \mu[t-j[q]],} 
and analoguesly for the model with logarithmic link by
    \deqn{\nu_t = \mu_t + \eta_1 X_{t,1} + \ldots + \eta_r X_{t,r},}{
    \nu[t] = \mu[t] + \eta[1] X[t,1] + \ldots + \eta[r] X[t,r],}
    \deqn{\mu_t = \beta_0 + \beta_1 \log(Z_{t-i_1}+1) + \ldots + \beta_p \log(Z_{t-i_p}+1)
    + \alpha_1 \mu{t-j_1} + \ldots  + \alpha_q \mu{t-j_q}.}{
    \mu[t] = \beta[0] + \beta[1] \log(Z[t-i[1]]+1) + \ldots + \beta[p] \log(Z[t-i[p]+1)
    + \alpha[1] \mu[t-j[1]] + \ldots  + \alpha[q] \mu[t-j[q]].} 
This is described in more detail by Liboschik et al. (2014) for the case of deterministic covariates for modelling interventions.
It is also possible to model a combination of external and internal covariates, which can be defined straightforwardly by adding each covariate either to the linear predictor \eqn{\nu_t}{\nu[t]} itself (for an internal effect) or to \eqn{\mu_t}{\mu[t]} defined above (for an external effect).

}

\value{
An object of class \code{"tsglm"}, which is a list with at least the following elements:

  \item{coefficients}{a named vector of the maximum likelihood estimated coefficients, which can be extracted by the \code{\link{coef}} method.}
  
  \item{start}{a named vector of the start estimation for the coefficients.}
  
  \item{residuals}{a vector of residuals, which can be extracted by the \code{\link[=residuals.tsglm]{residuals}} method.}

  \item{fitted.values}{the fitted values, which can be extracted by the \code{\link{fitted}} method.}

  \item{linear.predictors}{the linear fit on link scale.}
  
  \item{response}{a vector of the response values (this is usually the original time series but  possibly without the first few observations used for initialization if argument \code{init.drop=TRUE}).}

  \item{logLik}{the log-likelihood of the fitted model, which can be extracted by the \code{\link{logLik}} method. This is the complete log-likelihood including all constant terms. It is based on \code{n_eff} observations (see below).}

  \item{score}{the score vector at the maximum likelihood estimation.}

  \item{info.matrix}{the information matrix at the maximum likelihood estimation assuming a Poisson distribution.}

  \item{info.matrix_corrected}{the information matrix at the maximum likelihood estimation assuming the distribution specified in \code{distr}.}

  \item{call}{the matched call.}

  \item{n_obs}{the number of observations.}
  
  \item{n_eff}{the effective number of observations used for maximum likelihood estimation (might be lower than \code{n_obs} if argument \code{init.drop=TRUE}).}

  \item{ts}{the original time series.}

  \item{model}{the model specification.}
  
  \item{xreg}{the given covariates.}
  
  \item{distr}{a character giving the fitted conditional distribution.}
  
  \item{distrcoefs}{a named vector of the estimated additional coefficients specifying the conditional distribution. Is \code{NULL} in case of a Poisson distribution.}
  
  \item{sigmasq}{the estimated overdispersion coefficient. Is zero in case of a Poisson distribution.}

The function \code{tsglm.meanfit} has the same output except the elements \code{distr}, \code{distrcoefs} and \code{sigmasq}. In addition, they return the following list elements:

  \item{inter}{some details on the intermediate estimation of the coefficients as returned by \code{\link{constrOptim}} or \code{\link{optim}}.} 

  \item{final}{some details on the final estimation of the coefficients as returned by \code{\link{constrOptim}} or \code{\link{optim}}.}

  \item{durations}{named vector of the durations of the model fit (in seconds).}
  
  \item{outerscoreprod}{array of outer products of score vectors at each time point.} 

}

\author{
Tobias Liboschik, Philipp Probst, Konstantinos Fokianos and Roland Fried
}

\references{
Christou, V. and Fokianos, K. (2014) Quasi-likelihood inference for negative binomial time series models. \emph{Journal of Time Series Analysis} \bold{35(1)}, 55--78, \url{http://dx.doi.org/10.1002/jtsa.12050}.

Christou, V. and Fokianos, K. (2015) Estimation and testing linearity for non-linear mixed poisson autoregressions. \emph{Electronic Journal of Statistics} \bold{9}, 1357--1377, \url{http://dx.doi.org/10.1214/15-EJS1044}.

Ferland, R., Latour, A. and Oraichi, D. (2006) Integer-valued GARCH process. \emph{Journal of Time Series Analysis} \bold{27(6)}, 923--942, \url{http://dx.doi.org/10.1111/j.1467-9892.2006.00496.x}.

Fokianos, K. and Fried, R. (2010) Interventions in INGARCH processes. \emph{Journal of Time Series Analysis} \bold{31(3)}, 210--225, \url{http://dx.doi.org/10.1111/j.1467-9892.2010.00657.x}.

Fokianos, K., and Fried, R. (2012) Interventions in log-linear Poisson autoregression. \emph{Statistical Modelling} \bold{12(4)}, 299--322. \url{http://dx.doi.org/10.1177/1471082X1201200401}.

Fokianos, K., Rahbek, A. and Tjostheim, D. (2009) Poisson autoregression. \emph{Journal of the American Statistical Association} \bold{104(488)}, 1430--1439, \url{http://dx.doi.org/10.1198/jasa.2009.tm08270}.

Fokianos, K. and Tjostheim, D. (2011) Log-linear Poisson autoregression. \emph{Journal of Multivariate Analysis} \bold{102(3)}, 563--578, \url{http://dx.doi.org/10.1016/j.jmva.2010.11.002}.

Liboschik, T., Fokianos, K. and Fried, R. (2017) tscount: An R package for analysis of count time series following generalized linear models. \emph{Journal of Statistical Software} \bold{82(5)}, 1--51, \url{http://dx.doi.org/10.18637/jss.v082.i05}.
}

\seealso{
S3 methods \code{\link{print}}, \code{\link[=summary.tsglm]{summary}}, \code{\link[=residuals.tsglm]{residuals}}, \code{\link[=plot.tsglm]{plot}}, \code{\link{fitted}}, \code{\link{coef}}, \code{\link[=predict.tsglm]{predict}}, \code{\link{logLik}}, \code{\link{vcov}}, \code{\link{AIC}}, \code{\link{BIC}} and \code{\link[=QIC.tsglm]{QIC}} for the class \code{"tsglm"}.
The S3 method \code{\link[=se.tsglm]{se}} computes the standard errors of the parameter estimates.
Additionally, there are the S3 methods \code{\link[=pit.tsglm]{pit}}, \code{\link[=marcal.tsglm]{marcal}} and \code{\link[=scoring.tsglm]{scoring}} for predictive model assessment.


S3 methods \code{\link[=interv_test.tsglm]{interv_test}}, \code{\link[=interv_detect.tsglm]{interv_detect}} and \code{\link[=interv_multiple.tsglm]{interv_multiple}} for tests and detection procedures for intervention effects.
\code{\link{tsglm.sim}} for simulation from GLM-type model for time series of counts. \code{\link{ingarch.mean}}, \code{\link{ingarch.var}} and \code{\link{ingarch.acf}} for calculation of analytical mean, variance and autocorrelation function of an INGARCH model (i.e. with identity link) without covariates.

Example time series of counts are \code{\link{campy}}, \code{\link{ecoli}}, \code{\link{ehec}}, \code{\link{influenza}}, \code{\link{measles}} in this package, \code{\link[gamlss.data]{polio}} in package \code{gamlss.data}.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
interventions <- interv_covariate(n=length(campy), tau=c(84, 100),
  delta=c(1, 0)) #detected by Fokianos and Fried (2010, 2012)
#Linear link function with Negative Binomial distribution:
campyfit <- tsglm(campy, model=list(past_obs=1, past_mean=13),
  xreg=interventions, distr="nbinom")
campyfit
plot(campyfit)

###Road casualties in Great Britain (see help("Seatbelts"))
timeseries <- Seatbelts[, "VanKilled"]
regressors <- cbind(PetrolPrice=Seatbelts[, c("PetrolPrice")],
                    linearTrend=seq(along=timeseries)/12)
#Logarithmic link function with Poisson distribution:                    
seatbeltsfit <- tsglm(ts=timeseries, link="log",
  model=list(past_obs=c(1, 12)), xreg=regressors, distr="poisson")
summary(seatbeltsfit)
}
