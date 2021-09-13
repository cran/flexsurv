##' Simulate datasets from a fitted flexsurvreg model
##'
##' @param object Object returned by flexsurvreg
##'
##' @param nsim Number of datasets to simulate
##'
##' @param seed Random number seed. This is returned with the result of this
##'   function as described in \code{\link{simulate}}.
##'
##' @param censtime Right-censoring time, or vector of right-censoring times of
##'   the same size as the data used to fit the model.
##'   
##' @param vectorised Set to \code{FALSE} if the fitted model uses distribution
##'   functions from a package other than \pkg{flexsurv} and those functions are
##'   not vectorised.  By default, this function assumes that they are
##'   vectorised. Incomprehensible warnings may be printed and the simulation is
##'   not guaranteed to work if this assumption is wrong.  Vectorising will 
##'   generally make the simulation much faster. 
##'
##' @param ... Currently unused.
##'
##' @return A data frame with \code{nsim} pairs of columns named
##'   \code{"sim_1","event_1"} and so on, containing the simulated event or
##'   censoring times, and an indicator for whether the event was observed.
##'   
simulate.flexsurvreg <- function(object, nsim=1, seed=NULL, 
                                 censtime=NULL, vectorised=TRUE, ...){
  if (!is.null(seed)) set.seed(seed)
  dat <- model.frame(object)
  nobs <- nrow(dat)
  pars <- object$res[object$dlist$pars,"est"]
  pars <- t(as.matrix(pars))[rep(1, nobs),]
  beta <- t(object$res[object$covpars,"est",drop=FALSE])
  X <- model.matrix(object)
  if (object$ncovs > 0)
    pars <- add.covs(object, pars, beta, X)
  args <- c(as.list(as.data.frame(pars)), list(n=nobs))
  res <- as.data.frame(matrix(nrow=args$n, ncol=nsim))
  if (!is.null(seed)) attr(seed, "kind") <- as.list(RNGkind())
  seed_attr <- if (is.null(seed)) .Random.seed else seed
  if (!is.numeric(nsim) || (nsim<1)) 
    stop("`nsim` must be 1 or more")
  for (i in seq_len(nsim)) {
    if (vectorised)
      res[[i]] <- do.call(object$dfns$r, args)
    else {
      for (j in 1:nobs){
        argsj <- c(as.list(as.data.frame(pars)[j,]), list(n=1))
        res[[i]][j] <- do.call(object$dfns$r, argsj)
      }
    }
  }
  names(res) <- paste("time", seq_len(nsim), sep="_")
  if (is.null(censtime)) censtime <- Inf
  censtime <- rep(censtime, length.out = nobs)
  event <- as.data.frame(matrix(nrow=args$n, ncol=nsim))
  names(event) <- paste("event", seq_len(nsim), sep="_")
  
  for (i in seq_len(nsim)) {
    event[[i]] <- as.numeric(res[[i]] < censtime)
    res[[i]] <- ifelse(event[[i]], res[[i]], censtime)
  }
  col_order <- paste(rep(c("time","event"), nsim), 
                     rep(1:nsim, each=2), sep="_")
  res <- cbind(res, event)[,col_order]
  attr(res, "seed") <- seed_attr
  res
  
}

model.frame.flexsurvmix <- function(formula, par="time", k=1, ...) {
    x <- formula 
    if (par=="prob")
        res <- x$data$mf[[1]]
    else if (par=="time") {
        if (k %in% x$evnames)
            k <- match(k, x$evnames)
        else if (!(k %in% seq_len(x$K)))
            stop(sprintf("`k` should either be one of the event names or an integer in 1,...%s", x$K))
        res <- x$data$mf[[k]]
    }
    else stop("`par` should be either `prob` or `time`")
    res
}

model.matrix.flexsurvmix <- function(object, par="time", k=1, i=1, ...) {
  mf <- model.frame.flexsurvmix(object, par=par, k=k)
  if (par=="prob") {
    form <- object$pformula
  }
  else if (par=="time") 
    form <- object$all.formulae[[k]][[i]]
  model.matrix(form, mf)
}
