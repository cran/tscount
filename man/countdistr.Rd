\name{countdistr}
\alias{countdistr}
\alias{ddistr}
\alias{pdistr}
\alias{qdistr}
\alias{rdistr}
\alias{sddistr}
\alias{ardistr}
\alias{checkdistr}
\title{
Count Data Distributions 
}
\description{
Density, distribution function, quantile function, random generation, standard deviation and Anscombe residuals for some count data distributions. These auxiliary functions are used by several functions of the \code{tscount} package.
}
\usage{
ddistr(x, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, ...)
pdistr(q, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, ...)
qdistr(p, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, ...)
rdistr(n, meanvalue, distr=c("poisson", "nbinom"), distrcoefs)
sddistr(meanvalue, distr=c("poisson", "nbinom"), distrcoefs)
ardistr(response, meanvalue, distr=c("poisson", "nbinom"), distrcoefs)
checkdistr(distr=c("poisson", "nbinom"), distrcoefs)
}

\arguments{
  \item{x}{
   	vector of (non-negative integer) quantiles.
}
  \item{q}{
    vector of quantiles.
}
  \item{p}{
    vector of probabilities.
}
  \item{n}{
    positive integer value giving the number of random values to return.
}
  \item{response}{
    vector of true observations for calculation of residuals.
}
  \item{meanvalue}{
    non-negative numeric vector of means.
}
  \item{distr}{
    character value giving the distribution. Possible values are currently \code{"poisson"} (the default) for the \code{\link{Poisson}} distribution and \code{"nbinom"} for the \code{\link{NegBinomial}} distribution.
}
  \item{distrcoefs}{
    vector of additional distribution coefficients. For the Poisson distribution this argument can be omitted. For the negative binomial distribution it needs to be a vector of length one giving the value for the parameter \code{size} as specified in \code{\link{NegBinomial}}.
}
  \item{...}{
  additional arguments \code{logt}, \code{lower.tail} or \code{log.p} to be passed to other functions (see \code{\link{Poisson}} or \code{\link{NegBinomial}}).
  }
}

\details{
Basically, these function are wrappers for specific functions for the respective distribution. The function \code{ddistr} gives the density of the specified distribution, \code{pdistr} the distribution function, \code{qdistr} the quantile function and \code{rdistr} generates random deviates from this distribution. These functions are a generalisation of the respective functions where \code{distr} is replaced by either \code{pois} or \code{nbinom}. The function \code{sddistr} returns the standard deviation of the specified distribution. The function \code{ardistr} calculates Anscombe residuals for given values of the response. The function \code{checkdistr} is for verification of the arguments \code{distr} and \code{distrcoefs}. 
}

\author{
Tobias Liboschik
}


\seealso{
\code{\link{Poisson}} for the Poisson distribution and \code{\link{NegBinomial}} for the negative binomial distribution.

\code{\link{tsglm}} for fitting a more genereal GLM for time series of counts.
}
