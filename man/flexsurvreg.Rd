\encoding{latin1}
\name{flexsurvreg}
\alias{flexsurvreg}
\title{
  Flexible parametric regression for time-to-event data
}
\description{
  Parametric regression for time-to-event data using the generalized F
  and other flexible distributions.  Users may easily extend this
  function with their own survival distributions.
}
\usage{
flexsurvreg(formula, data, weights, subset, na.action, dist, inits,
            fixedpars=NULL, cl=0.95,...)
}
\arguments{
  \item{formula}{
    A formula expression in conventional R linear modelling
    syntax. The response must be a survival object as returned by the
    \code{\link{Surv}} function, and any covariates are given on the
    right-hand side.  For example,

    \code{Surv(time, dead) ~ age + sex}

    Only \code{Surv} objects
    of \code{type="right"} or \code{type="counting"}, corresponding to
    right-censored and/or left-truncated observations, are supported.

    If there are no covariates, specify \code{1} on the right hand side,
    for example \code{Surv(time, dead) ~ 1}.

    By default, covariates are placed on the ``location'' parameter of
    the distribution, typically the "scale" or "rate" parameter, through
    a linear model, or a log-linear model if this parameter must be
    positive.  This gives an accelerated failure time model or a
    proportional hazards model, depending on the distribution.

    Covariates can be placed on other parameters by using the name of
    the parameter as a function in the formula.  For example, in a
    Weibull model, the following expresses the scale parameter in terms
    of age and a treatment variable \code{treat}, and the shape
    parameter in terms of sex and treatment.

    \code{Surv(time, dead) ~ age + treat + shape(sex) + shape(treatment)}

  }
  \item{data}{
    A data frame in which to find variables supplied in \code{formula}.
    If not given, the variables should be in the working environment.
  }
  \item{weights}{
    Optional vector of case weights.
  }
  \item{subset}{
    Vector of integer or logicals specifying the subset of the
    observations to be used in the fit.
  }
  \item{na.action}{a missing-data filter function, applied after any
    'subset' argument has been used. Default is \code{options()$na.action}.
  }
  \item{dist}{ Either one of the following strings identifying a
    built-in distribution:
    \tabular{ll}{
      \code{"gengamma"} \tab Generalized gamma (stable parameterisation) \cr
      \code{"gengamma.orig"} \tab Generalized gamma (original parameterisation) \cr
      \code{"genf"} \tab Generalized F (stable parameterisation) \cr
      \code{"genf.orig"} \tab Generalized F (original parameterisation) \cr
      \code{"weibull"} \tab Weibull \cr
      \code{"gamma"} \tab Gamma \cr
      \code{"exp"} \tab Exponential \cr
      \code{"lnorm"} \tab Log-normal \cr
      \code{"gompertz"} \tab Gompertz \cr
    }

    or a list specifying a custom distribution.  See section ``Custom
    distributions'' below for how to construct this list.

    The parameterisations of the built-in distributions used here are
    the same as in their built-in distribution functions:
    \code{\link{dgengamma}}, \code{\link{dgengamma.orig}},
    \code{\link{dgenf}}, \code{\link{dgenf.orig}},
    \code{\link{dweibull}}, \code{\link{dgamma}}, \code{\link{dexp}},
    \code{\link{dlnorm}}, \code{\link{dgompertz}}, respectively.  The
    functions in base R are used where available, otherwise, they are
    provided in this package.

    Note that the Weibull parameterisation is different from that
    in \code{\link[survival]{survreg}}, instead it is
    consistent with \code{\link{dweibull}}.
    The \code{"scale"} reported by \code{\link[survival]{survreg}}
    is equivalent to \code{1/shape} as defined by
    \code{\link{dweibull}} and hence \code{\link{flexsurvreg}}.
    The first coefficient \code{(Intercept)} reported by
    \code{\link[survival]{survreg}} is equivalent to \code{log(scale)}
    in \code{\link{dweibull}} and \code{\link{flexsurvreg}}.

    Similarly in the exponential distribution, the rate, rather than the
    mean, is modelled on covariates.
  }
  \item{inits}{
    A numeric vector giving initial values for each unknown parameter.
    If not specified, default initial values are chosen from a simple
    summary of the uncensored survival time, for example the mean
    is often used to initialize scale parameters.  See the object
    \code{flexsurv.dists} in the source for the exact methods used.
    If the likelihood surface may be uneven, it is advised to run
    the optimisation starting from various different initial values
    to ensure convergence to the true global maximum.
  }
  \item{fixedpars}{Vector of indices of parameters whose values will be
    fixed at their initial values during the optimisation.  The indices
    are ordered with parameters of the baseline distribution coming
    first, followed by covariate effects.  For example, in a stable
    generalized Gamma model with two covariates, to fix the third
    of three generalized gamma parameters (the shape \code{Q},
    see the help for \code{\link{GenGamma}}) and the second covariate,
    specify \code{fixedpars = c(3, 5)}
  }
  \item{cl}{ Width of symmetric confidence intervals for maximum
    likelihood estimates, by default 0.95.
  }
  \item{...}{Optional arguments to the general-purpose \R
    optimisation routine \code{\link{optim}}.  For example, the BFGS
    optimisation algorithm is the default in \code{\link{flexsurvreg}},
    but this can be changed, for example to \code{method="Nelder-Mead"} which can be
    more robust to poor initial values.
    If the optimisation fails to converge, consider normalising the
    problem using, for example, \code{control=list(fnscale = 2500)}, for
    example, replacing 2500 by a number of the order of magnitude of the
    likelihood. If 'false' convergence is reported with a non-positive-definite Hessian,
    then consider tightening the tolerance criteria for convergence. If
    the optimisation takes a long time, intermediate steps can be
    printed using the \code{trace} argument of the control list. See
    \code{\link{optim}} for details.
  }
}
\details{
  Parameters are estimated by maximum likelihood using the
  algorithms available in the standard R \code{\link{optim}} function.
  Parameters defined to be positive are estimated on the log scale.
  Confidence intervals are estimated from the Hessian at the maximum,
  and transformed back to the original scale of the parameters.

  The usage of \code{\link{flexsurvreg}} is intended to be as similar as possible to
  \code{\link[survival]{survreg}} in the \pkg{survival} package.
}
\value{
  A list of class \code{"flexsurvreg"} with the following elements.

  \item{call}{A copy of the function call, for use in post-processing.}
\item{dlist}{List defining the survival distribution used.}
\item{res}{Matrix of maximum likelihood estimates and confidence
  limits, with parameters on their natural scales.}
\item{res.t}{Matrix of maximum likelihood estimates and confidence
  limits, with parameters all transformed to the real line.  The \code{\link{coef}},
  \code{\link{vcov}} and \code{\link{confint}} methods for
  \code{flexsurvreg} objects work on this scale.}
\item{loglik}{Log-likelihood}
\item{AIC}{Akaike's information criterion (-2*log likelihood + 2*number of
  estimated parameters)}
}
\references{
  Jackson, C. H. and Sharples, L. D. and Thompson, S. G. (2010)
  Survival models in health economic evaluations: balancing fit and
  parsimony to improve prediction. International Journal of Biostatistics
  6(1):Article 34.

  Cox, C. (2008) The generalized \eqn{F} distribution: An umbrella for
  parametric survival analysis.  Statistics in Medicine 27:4301-4312.

  Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007)
  Parametric survival analysis and taxonomy of hazard functions for the
  generalized gamma distribution.  Statistics in Medicine 26:4252-4374
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\section{Custom distributions}{
  \code{\link{flexsurvreg}} is intended to be easy to extend to handle
  new distributions.  To define a new distribution for use in
  \code{\link{flexsurvreg}}, construct a list with the following
  elements:

  \describe{
    \item{\code{name}:}{A string naming the distribution.
      If this is called \code{"dist"}, for example, then there must be a
      function called \code{ddist} in the working environment which
      defines the probability density,
      and a function called \code{pdist} which defines the probability
      distribution or cumulative density.  These functions may be in an add-on
      package (see below for an example) or may be user-written.
      Arguments other than parameters must be named in the conventional
      way -- for example \code{x} for the first argument of the density
      function, as in \code{\link{dnorm}(x, ...)} and \code{q} for the first
      argument of the probability function.
    }
    \item{\code{pars}:}{Vector of strings naming the parameters of the
      distribution. These must be the same names as the arguments of the
      density and probability functions.
    }
    \item{\code{location}:}{Name of the parameter which can be modelled as a linear
      function of covariates, possibly after transformation.}
    \item{\code{transforms}:}{Vector of R functions which transform the range of
    values taken by each parameter onto the real line.  For example,
    \code{c(log, log)} for a distribution with two positive parameters.}
    \item{\code{inv.transforms}:}{Vector of R functions defining the corresponding
      inverse transformations.}
    \item{\code{inits}:}{A function of the uncensored survival times \code{t},
      which returns a vector of reasonable initial values for maximum
      likelihood estimation of each parameter.  For example,
      \code{function(t){ c(1, mean(t)) }} will always initialize the first
      of two parameters at 1, and the second (a scale
      parameter, for instance) at the mean of \code{t}.
    }
  }

  For example, suppose we want to use a log-logistic survival
  distribution. This is available in the CRAN package \pkg{eha}, which
  provides conventionally-defined density and probability functions called
  \code{\link[eha]{dllogis}} and \code{\link[eha]{pllogis}}.  See the
  Examples below for the custom list in this case, and the
  subsequent command to fit the model.
}
\seealso{
  \code{\link{flexsurvspline}} for flexible survival modelling using the
  spline model of Royston and Parmar.

  \code{\link{plot.flexsurvreg}} and \code{\link{lines.flexsurvreg}} to
  plot fitted survival, hazards and cumulative hazards from models fitted
  by \code{\link{flexsurvreg}} and \code{\link{flexsurvspline}}.
}
\examples{
data(ovarian)
## Compare generalized gamma fit with Weibull
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gengamma")
fitg
fitw <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="weibull")
fitw
plot(fitg)
lines(fitw, col="blue", lwd.ci=1, lty.ci=1)
## Identical AIC, probably not enough data in this simple example for a
## very flexible model to be worthwhile.

## Custom distribution
library(eha)  ## make "dllogis" and "pllogis" available to the working environment
custom.llogis <- list(name="llogis",
                      pars=c("shape","scale"),
                      location="scale",
                      transforms=c(log, log),
                      inv.transforms=c(exp, exp),
                      inits=function(t){ c(1, median(t)) })
fitl <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.llogis)
fitl
lines(fitl, col.fit="purple", col.ci="purple")
}
\keyword{models}
\keyword{survival}
