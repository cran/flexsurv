% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spline.R
\name{flexsurvspline}
\alias{flexsurvspline}
\title{Flexible survival regression using the Royston/Parmar spline model.}
\usage{
flexsurvspline(
  formula,
  data,
  weights,
  bhazard,
  rtrunc,
  subset,
  k = 0,
  knots = NULL,
  bknots = NULL,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  ...
)
}
\arguments{
\item{formula}{A formula expression in conventional R linear modelling
syntax. The response must be a survival object as returned by the
\code{\link[survival]{Surv}} function, and any covariates are given on the right-hand
side.  For example,

\code{Surv(time, dead) ~ age + sex}

specifies a model where the log cumulative hazard (by default, see
\code{scale}) is a linear function of the covariates \code{age} and
\code{sex}.

If there are no covariates, specify \code{1} on the right hand side, for
example \code{Surv(time, dead) ~ 1}.

Time-varying covariate effects can be specified using the method described
in \code{\link{flexsurvreg}} for placing covariates on ancillary
parameters.  The ancillary parameters here are named \code{gamma1},
\ldots{}, \code{gammar} where \code{r} is the number of knots \code{k} plus
one (the ``degrees of freedom'' as defined by Royston and Parmar).  So for
the default Weibull model, there is just one ancillary parameter
\code{gamma1}.

Therefore a model with one internal spline knot, where the equivalents of
the Weibull shape and scale parameters, but not the higher-order term
\code{gamma2}, vary with age and sex, can be specified as:

\code{Surv(time, dead) ~ age + sex + gamma1(age) + gamma1(sex)}

or alternatively (and more safely, see \code{flexsurvreg})

\code{Surv(time, dead) ~ age + sex, anc=list(gamma1=~age + sex)}

\code{Surv} objects of \code{type="right"},\code{"counting"},
\code{"interval1"} or \code{"interval2"} are supported, corresponding to
right-censored, left-truncated or interval-censored observations.}

\item{data}{A data frame in which to find variables supplied in
\code{formula}.  If not given, the variables should be in the working
environment.}

\item{weights}{Optional variable giving case weights.}

\item{bhazard}{Optional variable giving expected hazards for relative
survival models.}

\item{rtrunc}{Optional variable giving individual right-truncation times (see \code{\link{flexsurvreg}}).  Note that these models can suffer from weakly identifiable parameters and
badly-behaved likelihood functions, and it is advised to compare
convergence for different initial values by supplying different
\code{inits} arguments to \code{flexsurvspline}.}

\item{subset}{Vector of integers or logicals specifying the subset of the
observations to be used in the fit.}

\item{k}{Number of knots in the spline. The default \code{k=0} gives a
Weibull, log-logistic or lognormal model, if \code{"scale"} is
\code{"hazard"}, \code{"odds"} or \code{"normal"} respectively.  \code{k}
is equivalent to \code{df-1} in the notation of \code{stpm} for Stata.  The
knots are then chosen as equally-spaced quantiles of the log uncensored
survival times, for example, at the median with one knot, or at the 33\%
and 67\% quantiles of log time (or time, see \code{"timescale"}) with two
knots.  To override this default knot placement, specify \code{knots}
instead.}

\item{knots}{Locations of knots on the axis of log time (or time, see
\code{"timescale"}).  If not specified, knot locations are chosen as
described in \code{k} above.  Either \code{k} or \code{knots} must be
specified. If both are specified, \code{knots} overrides \code{k}.}

\item{bknots}{Locations of boundary knots, on the axis of log time (or
time, see \code{"timescale"}).  If not supplied, these are are chosen as
the minimum and maximum log death time.}

\item{scale}{If \code{"hazard"}, the log cumulative hazard is modelled as a
spline function.

If \code{"odds"}, the log cumulative odds is modelled as a spline function.

If \code{"normal"}, \eqn{-\Phi^{-1}(S(t))}{-InvPhi(S(t))} is modelled as a
spline function, where \eqn{\Phi^{-1}()}{InvPhi()} is the inverse normal
distribution function \code{\link{qnorm}}.}

\item{timescale}{If \code{"log"} (the default) the log cumulative hazard
(or alternative) is modelled as a spline function of log time.  If
\code{"identity"}, it is modelled as a spline function of time, however
this model would not satisfy the desirable property that the cumulative hazard
(or alternative) should approach 0 at time zero.}

\item{spline}{\code{"rp"} to use the natural cubic spline basis
  described in Royston and Parmar.

  \code{"splines2ns"} to use the alternative natural cubic spline
  basis from the \code{splines2} package (Wang and Yan 2021),
  which may be better behaved due to the basis being orthogonal.}

\item{...}{Any other arguments to be passed to or through
\code{\link{flexsurvreg}}, for example, \code{anc}, \code{inits},
\code{fixedpars}, \code{weights}, \code{subset}, \code{na.action}, and any
options to control optimisation.  See \code{\link{flexsurvreg}}.}
}
\value{
A list of class \code{"flexsurvreg"} with the same elements as
described in \code{\link{flexsurvreg}}, and including extra components
describing the spline model.  See in particular:

\item{k}{Number of knots.} \item{knots}{Location of knots on the log time
axis.} \item{scale}{The \code{scale} of the model, hazard, odds or normal.}
\item{res}{Matrix of maximum likelihood estimates and confidence limits.
Spline coefficients are labelled \code{"gamma..."}, and covariate effects
are labelled with the names of the covariates.

Coefficients \code{gamma1,gamma2,...} here are the equivalent of
\code{s0,s1,...} in Stata \code{streg}, and \code{gamma0} is the equivalent
of the \code{xb} constant term.  To reproduce results, use the
\code{noorthog} option in Stata, since no orthogonalisation is performed on
the spline basis here.

In the Weibull model, for example, \code{gamma0,gamma1} are
\code{-shape*log(scale), shape} respectively in \code{\link{dweibull}} or
\code{\link{flexsurvreg}} notation, or (\code{-Intercept/scale},
\code{1/scale}) in \code{\link[survival]{survreg}} notation.

In the log-logistic model with shape \code{a} and scale \code{b} (as in
\code{\link[eha:Loglogistic]{eha::dllogis}} from the \pkg{eha} package), \code{1/b^a} is
equivalent to \code{exp(gamma0)}, and \code{a} is equivalent to
\code{gamma1}.

In the log-normal model with log-scale mean \code{mu} and standard
deviation \code{sigma}, \code{-mu/sigma} is equivalent to \code{gamma0} and
\code{1/sigma} is equivalent to \code{gamma1}.  } \item{loglik}{The
maximised log-likelihood.  This will differ from Stata, where the sum of
the log uncensored survival times is added to the log-likelihood in
survival models, to remove dependency on the time scale.}
}
\description{
Flexible parametric modelling of time-to-event data using the spline model
of Royston and Parmar (2002).
}
\details{
This function works as a wrapper around \code{\link{flexsurvreg}} by
dynamically constructing a custom distribution using
\code{\link{dsurvspline}}, \code{\link{psurvspline}} and
\code{\link{unroll.function}}.

In the spline-based survival model of Royston and Parmar (2002), a
transformation \eqn{g(S(t,z))} of the survival function is modelled as a
natural cubic spline function of log time \eqn{x = \log(t)}{x = log(t)}
plus linear effects of covariates \eqn{z}.

\deqn{g(S(t,z)) = s(x, \bm{\gamma}) + \bm{\beta}^T \mathbf{z}}{g(S(t,z)) =
s(x, gamma) + beta^T z}

The proportional hazards model (\code{scale="hazard"}) defines
\eqn{g(S(t,\mathbf{z})) = \log(-\log(S(t,\mathbf{z}))) =
\log(H(t,\mathbf{z}))}{g(S(t,z)) = log(-log(S(t,z))) = log(H(t,z))}, the
log cumulative hazard.

The proportional odds model (\code{scale="odds"}) defines
\eqn{g(S(t,\mathbf{z})) }{g(S(t,z)) = log(1/S(t,z) - 1)}\eqn{ =
\log(S(t,\mathbf{z})^{-1} - 1)}{g(S(t,z)) = log(1/S(t,z) - 1)}, the log
cumulative odds.

The probit model (\code{scale="normal"}) defines \eqn{g(S(t,\mathbf{z})) =
}{g(S(t,z)) = -InvPhi(S(t,z))}\eqn{ -\Phi^{-1}(S(t,\mathbf{z}))}{g(S(t,z))
= -InvPhi(S(t,z))}, where \eqn{\Phi^{-1}()}{InvPhi()} is the inverse normal
distribution function \code{\link{qnorm}}.

With no knots, the spline reduces to a linear function, and these models
are equivalent to Weibull, log-logistic and lognormal models respectively.

The spline coefficients \eqn{\gamma_j: j=1, 2 \ldots }{gamma_j: j=1, 2
\ldots}, which are called the "ancillary parameters" above, may also be
modelled as linear functions of covariates \eqn{\mathbf{z}}, as

\deqn{\gamma_j(\mathbf{z}) = \gamma_{j0} + \gamma_{j1}z_1 + \gamma_{j2}z_2
+ ... }{gamma_j(z) = \gamma_{j0} + \gamma_{j1}z_1 + \gamma_{j2}z_2 + ... }

giving a model where the effects of covariates are arbitrarily flexible
functions of time: a non-proportional hazards or odds model.

Natural cubic splines are cubic splines constrained to be linear beyond
boundary knots \eqn{k_{min},k_{max}}{kmin,kmax}.  The spline function is
defined as

\deqn{s(x,\boldsymbol{\gamma}) = \gamma_0 + \gamma_1 x + \gamma_2 v_1(x) + \ldots +
}{s(x,gamma) = gamma0 + gamma1 x + gamma2 v1(x) + ... + gamma_{m+1}
vm(x)}\deqn{ \gamma_{m+1} v_m(x)}{s(x,gamma) = gamma0 + gamma1 x + gamma2
v1(x) + ... + gamma_{m+1} vm(x)}

where \eqn{v_j(x)}{vj(x)} is the \eqn{j}th basis function

\deqn{v_j(x) = (x - k_j)^3_+ - \lambda_j(x - k_{min})^3_+ - (1 - }{vj(x) =
(x - kj)^3_+ - \lambda_j(x - kmin)^3_+ - (1 -\lambda_j) (x -
kmax)^3_+}\deqn{ \lambda_j) (x - k_{max})^3_+}{vj(x) = (x - kj)^3_+ -
\lambda_j(x - kmin)^3_+ - (1 -\lambda_j) (x - kmax)^3_+}

\deqn{\lambda_j = \frac{k_{max} - k_j}{k_{max} - k_{min}}}{\lambda_j =
(kmax - kj) / (kmax - kmin)}

and \eqn{(x - a)_+ = max(0, x - a)}.
}
\examples{

## Best-fitting model to breast cancer data from Royston and Parmar (2002)
## One internal knot (2 df) and cumulative odds scale

spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="odds")

## Fitted survival

plot(spl, lwd=3, ci=FALSE)

## Simple Weibull model fits much less well

splw <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=0, scale="hazard")
lines(splw, col="blue", ci=FALSE)

## Alternative way of fitting the Weibull

\dontrun{
splw2 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")
}

}
\references{
Royston, P. and Parmar, M. (2002).  Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data, with application to prognostic modelling and estimation of treatment
effects. Statistics in Medicine 21(1):2175-2197.

Wang W, Yan J (2021). Shape-Restricted Regression Splines with R
Package splines2. Journal of Data Science, 19(3), 498-517.

Jackson, C. (2016). flexsurv: A Platform for Parametric Survival Modeling
in R. Journal of Statistical Software, 70(8), 1-33.
doi:10.18637/jss.v070.i08
}
\seealso{
\code{\link{flexsurvreg}} for flexible survival modelling using
general parametric distributions.

\code{\link{plot.flexsurvreg}} and \code{\link{lines.flexsurvreg}} to plot
fitted survival, hazards and cumulative hazards from models fitted by
\code{\link{flexsurvspline}} and \code{\link{flexsurvreg}}.
}
\author{
Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\keyword{models}
\keyword{survival}
