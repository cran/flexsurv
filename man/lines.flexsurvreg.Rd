\name{lines.flexsurvreg}
\title{Add fitted flexible survival curves to a plot}
\alias{lines.flexsurvreg}
\concept{Survival}
\description{
  Add fitted survival curves from a \code{\link{flexsurvreg}} model fit
  to an existing plot.

  Equivalent to \code{\link{plot.flexsurvreg}(...,add=TRUE)}.
}
\usage{
\method{lines}{flexsurvreg}(x, X=NULL, 
        col.fit="red", lty.fit=1, lwd.fit=2,...)
}
\arguments{

  \item{x}{Output from \code{\link{flexsurvreg}}, representing a fitted
    survival model object.}

  \item{X}{Matrix of covariate values to produce fitted survival curves
    for, as described in \code{\link{plot.flexsurvreg}}. }

  \item{col.fit}{ Colour of the fitted survival curve(s). }
  \item{lty.fit}{ Line type of the fitted survival curve(s). }
  \item{lwd.fit}{ Line width of the fitted survival curve(s). }

  \item{...}{Other arguments to be passed to the generic
    \code{\link{plot}} and \code{\link{lines}} functions.}
}
\seealso{
  \code{\link{flexsurvreg}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
\keyword{aplot}
