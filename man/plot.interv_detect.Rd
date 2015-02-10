\name{plot.interv_detect}
\alias{plot.interv_detect}
\title{
Plot Test Statistic of Intervention Detection Procedure for Count Time Series Following Generalised Linear Models
}
\description{
Provides a plot of the test statistics of a test on an intervention in GLM-type count time series (as returned by \code{interv_detect.tsglm}) against time.
}
\usage{
\method{plot}{interv_detect}(x, ...)
}
\arguments{
  \item{x}{
  an object of class \code{"interv_detect"}, usually a result of a call to \code{\link{interv_detect.tsglm}}.
}
  \item{...}{
  additional arguments to be passed to function \code{\link{plot}}.
}
}

\author{
Tobias Liboschik and Philipp Probst
}

\seealso{
\code{\link[=interv_detect.tsglm]{interv_detect}} for detecting an intervention effect in GLM-type count time series and \code{\link{tsglm}} for fitting such a model.
}

\examples{
\dontrun{
###Campylobacter infections in Canada (see help("campy"))
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
campyfit_intervdetect <- interv_detect(fit=campyfit, taus=80:120,
                        delta=1, external=FALSE)
#This example runs about 20 minutes on a single processing unit,
#of course depending on its speed.
plot(campyfit_intervdetect)}
}

\keyword{Intervention detection}
