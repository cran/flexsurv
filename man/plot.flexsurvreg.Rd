\name{plot.flexsurvreg}
\title{Plots of fitted flexible survival models}
\alias{plot.flexsurvreg}
\concept{Survival}
\description{
  Plot fitted survival, cumulative hazard or hazard from a parametric
  model against nonparametric estimates to diagnose goodness-of-fit. 
}
\usage{
\method{plot}{flexsurvreg}(x, X=NULL, type="survival", tmax=NULL,
        col.obs="black", lty.obs=1, lwd.obs=1, 
        col.fit="red", lty.fit=1, lwd.fit=2, add=FALSE,...)
}
\arguments{

  \item{x}{Output from \code{\link{flexsurvreg}}, representing a fitted
    survival model object.}

  \item{X}{ Matrix of covariate values to produce fitted survival curves
    for. Columns represent different covariates, and rows represent
    multiple curves.  For example \code{matrix(c(1,2),nrow=2)} if there
    is only one covariate in the model, and we want to plot curves for
    covariate values of 1 and 2.

    For ``factor'' (categorical) covariates, the 0/1 values of the
    contrasts representing factor levels should be used.  For example,
    for a covariate \code{agegroup} with levels \code{20-29, 30-39,
    40-49, 50-59}, and baseline level \code{20-29}, there are three
    contrasts.  To plot curves for groups \code{20-29} and \code{40-49},
    supply \code{X = rbind(c(0,0,0), c(0,1,0))}, since all contrasts are
    zero for the baseline level, and the second contrast is ``turned
    on'' for the third level \code{40-49}.

    If there are only factor covariates in the model, then Kaplan-Meier
    curves are plotted for all distinct groups, and by default, fitted
    curves are also plotted for these groups.   To plot Kaplan-Meier and
    fitted curves for only a subset of groups, use \code{plot(survfit())}
    followed by \code{lines.flexsurvreg()}.

    If there are any continuous covariates, then a single population
    Kaplan-Meier curve is drawn. By default, a single fitted curve
    is drawn with the covariates set to their mean values in the data - 
    for categorical covariates, the means of the 0/1 indicator
    variables are taken.     
  }
  \item{type}{ \code{"survival"} for survival, to be plotted against
    Kaplan-Meier estimates from \code{\link[survival]{plot.survfit}}.

    \code{"cumhaz"} for cumulative hazard, plotted against transformed
    Kaplan-Meier estimates from \code{\link[survival]{plot.survfit}}.

    \code{"hazard"} for hazard, to be plotted against smooth nonparametric
    estimates from \code{\link[muhaz]{muhaz}}.  The nonparametric
    estimates tend to be unstable, and these plots are intended just to
    roughly indicate the shape of the hazards through time. 
  }
  \item{tmax}{ Maximum time point on the x-axis. By default, this is the
    maximum time observed in the data. }
  \item{col.obs}{ Colour of the nonparametric curve. }
  \item{lty.obs}{ Line type of the nonparametric curve. }
  \item{lwd.obs}{ Line width of the nonparametric curve. }

  \item{col.fit}{ Colour of the fitted parametric curve(s). }
  \item{lty.fit}{ Line type of the fitted parametric curve(s). }
  \item{lwd.fit}{ Line width of the fitted parametric curve(s). }

  \item{add}{If \code{TRUE}, add lines to an existing plot,
    otherwise new axes are drawn.} 

  \item{...}{Other options to be passed to 
    \code{\link{plot.survfit}} or \code{\link[muhaz]{muhaz}}, for
    example, to control the smoothness of the nonparametric hazard
    estimates. }
}
\note{
  Some standard plot arguments such as \code{"xlim","xlab"} may not work. 
  This function was designed as a quick check of model fit.  Users wanting
  publication-quality graphs are advised to set up an empty plot with
  the desired axes first (e.g. with \code{plot(...,type="n",...)}), then
  use suitable \code{\link{lines}} functions to add lines.   
}
\seealso{
  \code{\link{flexsurvreg}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
\keyword{hplot}
