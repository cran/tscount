\name{campy}
\alias{campy}
\title{
Campylobacter Infections Time Series
}
\description{
Time series with the number of cases of campylobacter infections in the north of the province Quebec (Canada) in four week intervals from January 1990 to the end of October 2000. It has 13 observations per year and 140 observations in total.
Campylobacterosis is an acute bacterial infectious disease attacking the digestive system.
}

\usage{
campy
}

\format{
A time series of class \code{"ts"}.
}

\source{
Ferland, R., Latour, A. and Oraichi, D. (2006) Integer-valued GARCH process. \emph{Journal of Time Series Analysis} \bold{27(6)}, 923--942, \url{http://dx.doi.org/10.1111/j.1467-9892.2006.00496.x}.
}

\seealso{
\code{\link{ecoli}}, \code{\link{ehec}}, \code{\link{influenza}}, \code{\link{measles}} in this package, \code{\link[gamlss.data]{polio}} in package \code{gamlss.data}
}

\examples{
plot(campy)

#Fit the INGARCH model used in Ferland et al. (2006):
campyfit <- tsglm(ts=campy, model=list(past_obs=1, past_mean=c(7,13)))
summary(campyfit)
plot(campyfit)
#Note that these parameter estimations differ from those obtained by
#Ferland et al. (2006). This might be due to a different initialisation
#of pre-sample values and different optimisation algorithms (they use
#Microsoft Excel Solver Macro). 
}

\keyword{Data}
