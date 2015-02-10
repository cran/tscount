\name{residuals.tsglm}
\alias{residuals.tsglm}
\title{
Residuals of a Generalised Linear Model for Time Series of Counts
}
\description{
Returns the residuals of a fitted GLM-type model for time series of counts.
}
\usage{
\method{residuals}{tsglm}(object, type = c("response", "pearson", "anscombe"), ...)
}

\arguments{
  \item{object}{
an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}.
}
  \item{type}{
character value giving the type of residuals which should be returned. Choose \code{type="response"} for raw residuals, \code{type="pearson"} for Pearson residuals and \code{type="anscombe"} for Anscombe residuals.
}
  \item{...}{
further arguments are currently ignored. Only for compatibility with generic function.
}
}
\details{
Computes a vector with the respective residuals of the fit given in argument \code{object}.
}

\value{
Numerical vector of the residuals.
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
campyfit_resid <- residuals(campyfit, type="pearson")
plot(campyfit_resid)
acf(campyfit_resid)
}

\keyword{Model assessment}
