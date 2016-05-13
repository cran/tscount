\name{QIC}
\alias{QIC}
\alias{QIC.tsglm}
\title{
Quasi Information Criterion of a Generalised Linear Model for Time Series of Counts
}
\description{
The function computes the Quasi Information Criterion (QIC) of a generalised linear model for time series of counts.
}
\usage{
\method{QIC}{tsglm}(object, ...)
}
\arguments{
  \item{object}{
an object of class \code{"tsglm"}.
}
  \item{...}{
additional arguments passed to \code{tscount:::tsglm.loglik}. These can be the arguments \code{init.method} and \code{init.drop} which are explained on the help page of the function \code{\link{tsglm}}.
}
}
\details{
The quasi information criterion (QIC) has been proposed by Pan (2001) as alternative to Akaike's information criterion (AIC) which is properly adjusted for regression analysis based on the generalized estimating equations (GEE).

This function computes the QIC of a generalised linear model for time series of counts. In case of models with the Poisson distribution the QIC has approximately the same value as the AIC. However, in case of models with another distribution it can be a more adequate alternative to the AIC.
}

\references{
Pan, W. (2001) Akaike's Information Criterion in Generalized Estimating Equations. \emph{Biometrics} \bold{57}, 120--125, \url{http://dx.doi.org/10.1111/j.0006-341X.2001.00120.x}.
}

\author{
Tobias Liboschik
}

\seealso{
\code{\link{tsglm}} for fitting a GLM for time series of counts.

\code{\link{AIC}} and \code{\link{BIC}} for other information criteria.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)), distr="nbinom")
QIC(campyfit)
AIC(campyfit)
}

\keyword{Model assessment}
