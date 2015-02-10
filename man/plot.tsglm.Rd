\name{plot.tsglm}
\alias{plot.tsglm}
\title{
Diagnostic Plots for a Fitted GLM-type Model for Time Series of Counts
}
\description{
Produces several diagnostic plots to asses the fit of a GLM-type model for time series of counts.
}
\usage{
\method{plot}{tsglm}(x, ask = TRUE, ...)
}
\arguments{
  \item{x}{
an object of class \code{"tsglm"}. Usually the result of a call to \code{\link{tsglm}}.
}
  \item{ask}{
logical value. If \code{TRUE} (and the R session is interactive) the user is asked for input, before a new figure is drawn (see \code{\link{devAskNewPage}}).
}
  \item{...}{
further arguments are currently ignored. Only for compatibility with generic function.
}
}
\details{
Produces plots of the acf of the Pearson residuals, the Pearson residuals plotted against time, a cumulative periodogramm of the Pearson residuals, a probability integral transform (PIT) histogram (see function \code{\link{pit}}) and a marginal calibration plot (see function \code{\link{marcal}}). The cumulative periodogramm is plotted with the function \code{\link{cpgram}} from package \code{MASS} and is omitted with a warning if this package is not available.
}
\author{
Tobias Liboschik and Philipp Probst
}

\seealso{
\code{\link{tsglm}} for fitting a GLM for time series of counts.
}

\examples{
###Campylobacter infections in Canada (see help("campy"))
interventions <- interv_covariate(n=length(campy), tau=c(84, 100),
  delta=c(1, 0)) #detected by Fokianos and Fried (2010, 2012)
#Linear link function with Negative Binomial distribution:
campyfit <- tsglm(campy, model=list(past_obs=1, past_mean=13),
  xreg=interventions, dist="nbinom")
plot(campyfit)
}

\keyword{Model assessment}
