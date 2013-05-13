## Deriv of loglik wrt transformed parameters p
## loglik(p|x) = sum(log(f(p|xobs))) + sum(log(S(p|xcens))) - sum(log(S(p|xtrunc)))
## dloglik/dp  = sum (df/dp / f(p)) | xobs) + sum(dS/dp / S(p) | xcens) - sum(dS/dp / S(p) | xtrunc)
##             = sum(dlogf/dp | xobs) + sum(dlogS/dp | xcens) - sum(dlogS/dp | xtrunc)

Dminusloglik.flexsurv <- function(optpars, Y, X=0, dlist, inits, trunc, fixedpars=NULL) {
    pars <- inits
    npars <- length(pars)
    pars[setdiff(1:npars, fixedpars)] <- optpars
    nbpars <- length(dlist$pars)
    pars <- as.list(pars)
    ncovs <- length(pars) - length(dlist$pars)
    if (ncovs > 0) {
        beta <- unlist(pars[(nbpars+1):npars])
        pars[[dlist$location]] <- pars[[dlist$location]] + X %*% beta
    }
    else pars[[dlist$location]] <- rep(pars[[dlist$location]], length(Y[,"stop"]))
    ddfn <- paste("DLd",dlist$name,sep="")
    dead <- Y[,"status"]==1
    ddcall <- list(t=Y[dead,"stop"], X=X[dead,,drop=FALSE], ncovs=ncovs)
    dsfn <- paste("DLS",dlist$name,sep="")
    dsccall <- list(t=Y[!dead,"stop"], X=X[!dead,,drop=FALSE], ncovs=ncovs)
    dstcall <- list(t=Y[,"start"], X=X, ncovs=ncovs)
    for (i in 1:nbpars)
        ddcall[[names(pars)[i]]] <-
            dsccall[[names(pars)[i]]] <-
                dstcall[[names(pars)[i]]] <-
                    dlist$inv.transforms[[i]](pars[[i]])
    ddcall[[dlist$location]] <- ddcall[[dlist$location]][dead]
    dsccall[[dlist$location]] <- dsccall[[dlist$location]][!dead]
    dd <- do.call(ddfn, ddcall)
    dscens <- do.call(dsfn, dsccall)
    dstrunc <- do.call(dsfn, dstcall)
    res <- - ( colSums(dd) + colSums(dscens) - colSums(dstrunc) )
    ## currently wastefully calculates derivs for fixed pars then discards them
    res[setdiff(1:npars, fixedpars)]
}

## Derivatives of log density and log survival probability wrt
## baseline parameters and covariate effects for various distributions
## No easy derivatives available for other distributions.

## Exponential

DLdexp <- function(t, rate, X, ncovs){
    res <- matrix(nrow=length(t), ncol=1 + ncovs) # generic
    ts <- 1 - t*rate
    res[,1] <- ts
    for (i in seq_len(ncovs))
        res[,1+i] <- X[,i]*res[,1] # generic?
    res
}

DLSexp <- function(t, rate, X, ncovs){
    res <- matrix(nrow=length(t), ncol=1 + ncovs) # generic
    res[,1] <- -t*rate
    for (i in seq_len(ncovs))
        res[,1+i] <- X[,i]*res[,1]
    res
}

## Weibull

DLdweibull <- function(t, shape, scale, X, ncovs){
    res <- matrix(nrow=length(t), ncol=2 + ncovs)
    tss <- (t/scale)^shape
    res[,1] <- 1/shape + log(t/scale) - log(t/scale)*tss
    res[,2] <- -1 - (shape-1) + shape*tss
    for (i in seq_len(ncovs))
        res[,2+i] <- X[,i]*res[,2]
    res
}

DLSweibull <- function(t, shape, scale, X, ncovs){
    res <- matrix(nrow=length(t), ncol=2 + ncovs)
    tss <- (t/scale)^shape
    res[,1] <- ifelse(t==0, 0, -log(t/scale)*tss)
    res[,2] <- tss*shape
    for (i in seq_len(ncovs))
        res[,2+i] <- X[,i]*res[,2]
    res
}

## Gompertz

DLdgompertz <- function(t, shape, rate, X, ncovs){
    res <- matrix(nrow=length(t), ncol=2 + ncovs)
    rs <- rate/shape*exp(shape*t)
    res[,1] <- if (shape==0) 0 else t + rs*(1/shape - t) - rate/shape^2
    res[,2] <- if (shape==0) 1 - rate*t else 1 - rs + rate/shape
    for (i in seq_len(ncovs))
        res[,2+i] <- X[,i]*res[,2]
    res
}

DLSgompertz <- function(t, shape, rate, X, ncovs){
    res <- matrix(nrow=length(t), ncol=2 + ncovs)
    rs <- rate/shape*exp(shape*t)
    res[,1] <- if (shape==0) 0 else rs*(1/shape - t) - rate/shape^2
    res[,2] <- if (shape==0) -rate*t else - rs + rate/shape
    for (i in seq_len(ncovs))
        res[,2+i] <- X[,i]*res[,2]
    res
}


## Spline

Dminusloglik.stpm <- function(optpars, knots, Y, X=0, inits, fixedpars=NULL, scale="hazard"){
    pars <- inits
    npars <- length(pars)
    pars[setdiff(1:npars, fixedpars)] <- optpars
    nk <- length(knots)
    gamma <- pars[1:nk] # always at least two of these, intercept + weibull par
    if (npars > nk) {
        beta <- pars[(nk+1):npars]
    }
    else {beta <- 0; X <- matrix(0, nrow=nrow(Y))}
    dead <- Y[,"status"]==1
    dd <- DLdsurvspline(Y[dead,"stop"], gamma, beta, X[dead,,drop=FALSE], knots, scale)
    dscens <- DLSsurvspline(Y[!dead,"stop"], gamma, beta, X[!dead,,drop=FALSE], knots, scale)
    dstrunc <- DLSsurvspline(Y[,"start"], gamma, beta, X[,,drop=FALSE], knots, scale)
    res <- - ( colSums(dd) + colSums(dscens) - colSums(dstrunc) )
    res[setdiff(1:npars, fixedpars)]
}

DLdsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard"){
    res <- matrix(nrow=length(t), ncol=length(gamma)+length(beta))
    b <- basis(knots, log(t))
    db <- dbasis(knots, log(t))
    eta <- b %*% gamma + as.numeric(X %*% beta)
    ds <- db %*% gamma
    for (i in seq_along(gamma)){
        if (scale=="hazard")
            res[,i] <- db[,i] / ds + b[,i] * (1 - exp(eta))
        else if (scale=="odds"){
            eeta <- 1 - 2*exp(eta)/(1 + exp(eta))
            res[,i] <- db[,i] / ds + b[,i] * eeta
        }
    }
    for (i in seq_along(beta)){
        if (scale=="hazard")
            res[,length(gamma)+i] <- X[,i] * (1 - exp(eta))
        else if (scale=="odds")
            res[,length(gamma)+i] <- X[,i] * eeta
    }
    res
}

DLSsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard"){
    res <- matrix(nrow=length(t), ncol=length(gamma)+length(beta))
    b <- basis(knots, log(t))
    eta <- b %*% gamma + as.numeric(X %*% beta)
    for (i in seq_along(gamma)){
        if (scale=="hazard")
            res[,i] <- ifelse(t==0, 0, - b[,i] * exp(eta))
        else if (scale=="odds"){
            eeta <- exp(eta)/(1 + exp(eta))
            res[,i] <- ifelse(t==0, 0, - b[,i] * eeta)
        }
    }
    for (i in seq_along(beta)){
        if (scale=="hazard")
            res[,length(gamma)+i] <- ifelse(t==0, 0, - X[,i] * exp(eta))
        else if (scale=="odds")
            res[,length(gamma)+i] <- ifelse(t==0, 0, - X[,i] * eeta)
    }
    res
}
