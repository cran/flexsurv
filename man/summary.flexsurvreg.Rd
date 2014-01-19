\name{summary.flexsurvreg}
\title{Summaries of fitted flexible survival models}
\alias{summary.flexsurvreg}
\concept{Survival}
\description{
  Return fitted survival, cumulative hazard or hazard
  at a series of times from a fitted \code{\link{flexsurvreg}} or
  \code{\link{flexsurvspline}} model.
}
\usage{
\method{summary}{flexsurvreg}(object, X=NULL, type="survival", t=NULL,
         start=NULL, B=1000, cl=0.95,...)
}
\arguments{
  \item{object}{Output from \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}, representing a fitted
    survival model object.}

  \item{X}{ Matrix of covariate values to produce fitted values
    for. Columns represent different covariates, and rows represent
    multiple combinations of covariate values.
    For example \code{matrix(c(1,2),nrow=2)} if there
    is only one covariate in the model, and we want survival for
    covariate values of 1 and 2.

    For ``factor'' (categorical) covariates, the values of the contrasts
    representing factor levels (as returned by the
    \code{\link{contrasts}} function) should be used.  For example, for
    a covariate \code{agegroup} specified as an unordered factor with
    levels \code{20-29, 30-39, 40-49, 50-59}, and baseline level
    \code{20-29}, there are three contrasts.  To return summaries for
    groups \code{20-29} and \code{40-49}, supply
    \code{X = rbind(c(0,0,0), c(0,1,0))},
    since all contrasts are zero for the baseline level, and the second
    contrast is ``turned on'' for the third level \code{40-49}.

    If there are only factor covariates in the model, then all distinct
    groups are used by default.

    If there are any continuous covariates, then a single summary is
    provided.  By default, this is with all
    covariates set to their mean values in the data - for categorical
    covariates, the means of the 0/1 indicator variables are taken.
  }
  \item{type}{ \code{"survival"} for survival probabilities.

    \code{"cumhaz"} for cumulative hazards.

    \code{"hazard"} for hazards.
  }
  \item{t}{ Times to calculate fitted values for. By default, these are the
    sorted unique observation (including censoring) times in the
    data. If the corresponding left-truncation times \code{start} are not
    supplied, then they all default to 0.
  }

  \item{start}{ Left-truncation times, defaults to those corresponding
  to the default \code{t} in the data. }

  \item{B}{Number of simulations from the normal asymptotic distribution
  of the estimates used to calculate confidence intervals.  Decrease
    for greater speed at the expense of accuracy, or set
  \code{B=0} to turn off calculation of CIs.}

  \item{cl}{Width of symmetric confidence intervals, relative to 1.}

  \item{...}{Further arguments passed to or from other methods.}
}
\value{
  A list with one element for each unique covariate value (if there are
  only categorical covariates) or one element (if there are no
  covariates or any continuous covariates).  Each of these elements
  is a matrix with one row for each time in \code{t}, giving the
  estimated survival (or cumulative hazard, or hazard) and 95\%
  confidence limits.  These list elements are named with the covariate
  names and values which define them.

  If there are multiple summaries, an additional list component named
  \code{X} contains a matrix with the exact values of contrasts (dummy
  covariates) defining each summary.

  The \code{\link{plot.flexsurvreg}} function can be used to quickly
  plot these model-based summaries against empirical summaries such as
  Kaplan-Meier curves, to diagnose model fit.

  Confidence intervals for models fitted with \code{\link{flexsurvreg}}
  are obtained by random sampling from the asymptotic normal
  distribution of the maximum likelihood estimates (see, e.g. Mandel
  (2013)).  For models fitted
  with \code{\link{flexsurvreg}}, intervals for the hazard are obtained
  in this way, whereas intervals for the survival and cumulative hazard
  are obtained analytically as in Royston and Parmar (2002).
}
\seealso{
  \code{\link{flexsurvreg}}, \code{\link{flexsurvspline}}.
}
\references{
  Royston, P. and Parmar, M. (2002).  Flexible parametric
  proportional-hazards and proportional-odds models for censored survival
  data, with application to prognostic modelling and estimation of
  treatment effects. Statistics in Medicine 21(1):2175-2197.

  Mandel, M. (2013). "Simulation based confidence intervals for
  functions with complicated derivatives." The American Statistician (in
  press).
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
