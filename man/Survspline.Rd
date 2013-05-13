\name{Survspline}
\alias{dsurvspline}
\alias{psurvspline}
\title{Royston/Parmar spline survival distribution}
\description{
  Probability density and distribution function for the
  Royston/Parmar spline model.
}
\usage{
dsurvspline(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0)
psurvspline(q, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0)
}
\arguments{
  \item{x,q}{Vector of times.}
  \item{gamma}{Vector of parameters describing the baseline spline
    function, as described in \code{\link{flexsurvspline}}.}
  \item{beta}{Vector of covariate effects.}
  \item{X}{Matrix of covariate values.}
  \item{knots}{Locations of knots on the axis of log time, supplied in
    increasing order.  Unlike in \code{\link{flexsurvspline}}, these
    include the two boundary knots.  If there are
    no additional knots, the boundary locations are not used.  If there
    are one or more additional knots, the boundary knots should be at
    or beyond the minimum and maximum values of the log times. 
    In \code{\link{flexsurvspline}} these are exactly at the minimum
    and maximum values.
  }
  \item{scale}{\code{"hazard"}, \code{"odds"}, or \code{"normal"}, as
    described in \code{\link{flexsurvspline}}.  With the default of
    no knots in addition to the boundaries, this model reduces to the
    Weibull, log-logistic and log-normal respectively.}
  \item{offset}{An extra constant to add to the linear predictor \eqn{\eta}{eta}.}
}
\value{
\code{dsurvspline} gives the density and \code{psurvspline} gives the distribution
function, as described in \code{\link{flexsurvspline}}.
}
\references{
  Royston, P. and Parmar, M. (2002).  Flexible parametric
  proportional-hazards and proportional-odds models for censored survival
  data, with application to prognostic modelling and estimation of
  treatment effects. Statistics in Medicine 21(1):2175-2197. 
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\seealso{
   \code{\link{flexsurvspline}}.
}
\examples{
## reduces to the weibull
regscale <- 0.786; cf <- 1.82
a <- 1/regscale; b <- exp(cf)
dweibull(1, shape=a, scale=b)
dsurvspline(1, gamma=c(log(1 / b^a), a)) # should be the same

## reduces to the log-normal
meanlog <- 1.52; sdlog <- 1.11
dlnorm(1, meanlog, sdlog) 
dsurvspline(1, gamma = c(-meanlog/sdlog, 1/sdlog), scale="normal")
# should be the same
}
\keyword{ distribution }
