\name{normboot.flexsurvreg}
\alias{normboot.flexsurvreg}
\title{
  Simulate from the asymptotic normal distribution of parameter estimates.
}
\description{
  Produce a matrix of alternative parameter estimates under sampling
  uncertainty, at covariate values supplied by the user.
  Used by \code{\link{summary.flexsurvreg}} for obtaining
  confidence intervals around functions of parameters.
}
\usage{
normboot.flexsurvreg(x, B, newdata = NULL, X = NULL,
                     transform = FALSE, raw=FALSE)
}

\arguments{
  \item{x}{
    A fitted model from \code{\link{flexsurvreg}}.
  }
  \item{B}{
    Number of samples.
  }
  \item{newdata}{
    Data frame or list containing the covariate values to
    evaluate the parameters at.   If there are covariates in the model, at
    least one of \code{newdata} or \code{X} must be supplied,
    unless \code{raw=TRUE}.
  }
  \item{X}{
    Alternative (less convenient) format for covariate values: a
    matrix with one row, with one column for each covariate or factor
    contrast.  Formed from all the "model matrices", one for each
    named parameter of the distribution, with intercepts
    excluded, \code{cbind}ed together. 
  }
  \item{transform}{
   \code{TRUE} if the results should be transformed to the real-line
   scale, typically by log if the parameter is defined as positive.  The
   default \code{FALSE} returns parameters on the natural scale.
 }
 \item{raw}{
   Return samples of the baseline parameters and the covariate effects,
   rather than the default of adjusting the baseline parameters for
   covariates.}
}
\value{
  If \code{newdata} includes only one covariate combination, a matrix will be
  returned with \code{B} rows, and one column for each named parameter
  of the survival distribution.

  If more than one covariate combination is requested
  (e.g. \code{newdata} is a data frame with more than one row), then a
  list of matrices will be returned, one for each covariate combination.
}
\references{
  Mandel, M. (2013). "Simulation based confidence intervals for
  functions with complicated derivatives." The American Statistician (in
  press).
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\seealso{
\code{\link{summary.flexsurvreg}}
}
\examples{
    fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
    normboot.flexsurvreg(fite, B=10, newdata=list(age=50))
    normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
    normboot.flexsurvreg(fite, B=10, newdata=list(age=0))  ## closer to...
    fite$res
}
\keyword{ models }
