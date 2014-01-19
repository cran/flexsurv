
flexsurv.dists <- list(
                       genf = list(
                       pars=c("mu","sigma","Q","P"),
                       location="mu",
                       transforms=c(identity, log, identity, log),
                       inv.transforms=c(identity, exp, identity, exp),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt), 0, 1)
                       }
                       ),
                       genf.orig = list(
                       pars=c("mu","sigma","s1","s2"),
                       location="mu",
                       transforms=c(identity, log, log, log),
                       inv.transforms=c(identity, exp, exp, exp),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt), 1, 1)
                       }
                       ),
                       gengamma = list(
                       pars=c("mu","sigma","Q"),
                       location="mu",
                       transforms=c(identity, log, identity),
                       inv.transforms=c(identity, exp, identity),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt), 0)
                       }
                       ),
                       gengamma.orig = list(
                       pars=c("shape","scale","k"),
                       location="scale",
                       transforms=c(log, log, log),
                       inv.transforms=c(exp, exp, exp),
                       inits=function(t){c(1, mean(t), 1)}
                       ),
                       exp = list(
                       pars=c("rate"),
                       location="rate",
                       transforms=c(log),
                       inv.transforms=c(exp),
                       inits=function(t){1 / mean(t)},
                       deriv=TRUE
                       ),
                       weibull = list(
                       pars=c("shape","scale"),
                       location="scale",
                       transforms=c(log, log),
                       inv.transforms=c(exp, exp),
                       inits = function(t){
                           lt <- log(t[t>0])
                           c(1, exp(mean(lt) + 0.572))
                       },
                       deriv=TRUE
                       ),
                       lnorm = list(
                       pars=c("meanlog","sdlog"),
                       location="meanlog",
                       transforms=c(identity, log),
                       inv.transforms=c(identity, exp),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt))
                       }
                       ),
                       gamma = list(
                       pars=c("shape","rate"),
                       location="rate",
                       transforms=c(log, log),
                       inv.transforms=c(exp, exp),
                       inits=function(t){
                           m=mean(t); v=var(t);
                           c(m^2/v, m/v)
                       }
                       ),
                       gompertz = list(
                       pars=c("shape","rate"),
                       location="rate",
                       transforms=c(identity, log),
                       inv.transforms=c(identity, exp),
                       inits=function(t){c(0.001,1 / mean(t))},
                       deriv=TRUE
                       )
                       )

minusloglik.flexsurv <- function(optpars, Y, X=0, weights, dlist, inits, mx, fixedpars=NULL) {
    pars <- inits
    npars <- length(pars)
    pars[setdiff(1:npars, fixedpars)] <- optpars
    nbpars <- length(dlist$pars)
    pars <- as.list(pars)
    if (npars > nbpars) {
        beta <- unlist(pars[(nbpars+1):npars])
        for (i in dlist$pars)
            pars[[i]] <- pars[[i]] + X[,mx[[i]],drop=FALSE] %*% beta[mx[[i]]]
#        pars[[dlist$location]] <- pars[[dlist$location]] + X %*% beta
    }
    pcall <- list(q=Y[,"stop"])
    dcall <- list(x=Y[,"stop"])
    tcall <- list(q=Y[,"start"])
    for (i in 1:nbpars){
        pcall[[names(pars)[i]]] <-
            dcall[[names(pars)[i]]] <-
                tcall[[names(pars)[i]]] <-
                    dlist$inv.transforms[[i]](pars[[i]])
#        if (is.infinite(dlist$inv.transforms[[i]](pars[[i]])))
#            return(Inf)
    }
    dcall$log <- TRUE
    probfn <- paste("p",dlist$name,sep="")
    densfn <- paste("d",dlist$name,sep="")
    ## Generic survival model likelihood
    dead <- Y[,"status"]==1
    logdens <- (do.call(densfn, dcall)*weights)[dead]
    prob <- (do.call(probfn, pcall))[!dead]
    pobs <- 1 - do.call(probfn, tcall) # prob of being observed = 1 unless left-truncated
    - ( sum(logdens) + sum(log(1 - prob)*weights[!dead]) - sum(log(pobs)*weights))
}

check.dlist <- function(dlist){
## TESTME
    if (is.null(dlist$name)) stop("\"name\" element of custom distribution list not given")
    if (!is.character(dlist$name)) stop("\"name\" element of custom distribution list should be a string")
    if (is.null(dlist$pars)) stop("Parameter names \"pars\" not given in custom distribution list")
    if (!is.character(dlist$pars)) stop("Parameter names \"pars\" should be a character vector")
    npars <- length(dlist$pars)
    if (is.null(dlist$location)) stop("Location parameter not given in custom distribution list")
    if (!(dlist$location %in% dlist$pars)) stop("Location parameter \"",dlist$location,"\" not in list of parameters")
    if (is.null(dlist$transforms)) stop("Transforms not given in custom distribution list")
    if (length(dlist$transforms) != npars) stop("Transforms vector of length ",length(dlist$transforms),", parameter names of length ",npars)
    if (is.null(dlist$inv.transforms)) stop("Inverse transformations not given in custom distribution list")
    if (length(dlist$inv.transforms) != npars) stop("Inverse transforms vector of length ",length(dlist$transforms),", parameter names of length ",npars)
    for (i in 1:npars){
        if (is.character(dlist$transforms[[i]])) dlist$transforms[[i]] <- get(dlist$transforms[[i]])
        if (is.character(dlist$inv.transforms[[i]])) dlist$inv.transforms[[i]] <- get(dlist$inv.transforms[[i]])
        if (!is.function(dlist$transforms[[i]])) stop("Transformation function for parameter ", i, " not defined")
        if (!is.function(dlist$inv.transforms[[i]])) stop("Inverse transformation function for parameter ", i, " not defined")
    }
    if (!is.function(dlist$inits)) stop("\"inits\" element of custom distribution list must be a function")
    res <- dlist$inits(1:10)
    if (!is.numeric(res) || (length(res) != npars))
        stop("\"inits\" function must return a numeric vector of length ", npars, " = number of parameters")
    dlist
}

flexsurvreg <- function(formula, data, weights, subset, na.action, dist, inits, fixedpars=NULL, cl=0.95, ...)
{
    call <- match.call()
    if (missing(dist)) stop("Distribution \"dist\" not specified")
    if (is.character(dist)) {
        match.arg(dist, names(flexsurv.dists))
        dlist <- flexsurv.dists[[dist]]
        dlist$name <- dist
    }
    else if (is.list(dist)) {
        dlist <- check.dlist(dist)
    }
    else stop("\"dist\" should be a string for a built-in distribution, or a list for a custom distribution")
    parnames <- dlist$pars
    ancnames <- setdiff(parnames, dlist$location)

    indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(call), nomatch = 0)
    if (indx[1] == 0)
        stop("A \"formula\" argument is required")
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    ## local environment to facilitate formulae for covariates on ancillary parameters.  Thanks to Milan Bouchet-Valat.
    tmpenv <- new.env(parent=environment(formula))
    f2 <- formula
    environment(f2) <- tmpenv
    temp[["formula"]] <- f2
    for (i in ancnames)
        assign(i, identity, envir=tmpenv)
    m <- eval(temp, parent.frame())
    Terms <- terms(formula, ancnames)
    X <- model.matrix(Terms, m)

    inds <- lapply(ancnames, function(nam) survival::untangle.specials(Terms, nam, 1:2)$terms)
    ass <- attributes(X)$assign
    ancidx <- lapply(inds, function(x){which(ass %in% x) - 1})

    names(ancidx) <- ancnames
    X <- X[,-1,drop=FALSE]
    locidx <- setdiff(seq_len(ncol(X)), unlist(ancidx))
    mx <- c(list(locidx), ancidx); names(mx)[1] <- dlist$location
    mx <- mx[parnames] # sort in original order

    weights <- model.extract(m, "weights")
    if (is.null(weights)) weights <- rep(1, nrow(X))
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    if (!(attr(Y, "type")  %in% c("right","counting")))
        stop("Survival object type \"", attr(Y, "type"), "\"", " not supported")
    if (attr(Y, "type") == "counting")
        Y <- cbind(Y, time=Y[,"stop"] - Y[,"start"]) # converts Y from Surv object to numeric matrix
    else Y <- cbind(Y, start=0, stop=Y[,"time"])

    dat <- list(Y=Y, X=X, Xraw=m[,-1,drop=FALSE])
    ncovs <- ncol(dat$Xraw)
    ncoveffs <- ncol(X)
    nbpars <- length(parnames) # number of baseline parameters
    npars <- nbpars + ncoveffs
    if (!missing(inits) && (!is.numeric(inits) || (length(inits) != npars)))
        stop("inits must be a numeric vector of length ",npars)
    if (missing(inits) || any(is.na(inits)))
        default.inits <- c(dlist$inits(Y[,"time"]*weights*length(Y[,"time"])/sum(weights)), rep(0,ncoveffs))
    if (missing(inits)) inits <- default.inits
    else if (any(is.na(inits))) inits[is.na(inits)] <- default.inits[is.na(inits)]
    for (i in 1:nbpars)
        inits[i] <- dlist$transforms[[i]](inits[i])
    outofrange <- which(is.nan(inits) | is.infinite(inits))
    if (any(outofrange)){
        plural <- if(length(outofrange) > 1) "s" else ""
        stop("Initial value", plural, " for parameter", plural, " ",
             paste(outofrange,collapse=","), " out of range")
    }
    cnames <- if(ncoveffs==0) NULL else colnames(X)
    names(inits) <- c(parnames, cnames)
    if (!is.null(fixedpars) && !is.logical(fixedpars) &&
        (!is.numeric(fixedpars) || any(!(fixedpars %in% 1:npars)))){
        dots <- if(npars>2) "...," else ""
        stop("fixedpars must be TRUE/FALSE or a vector of indices in 1,",dots,npars)
    }

    if ((is.logical(fixedpars) && fixedpars==TRUE) ||
        (is.numeric(fixedpars) && all(fixedpars == 1:npars))) {
        minusloglik <- minusloglik.flexsurv(inits, Y=Y, X=X, weights=weights, dlist=dlist, inits=inits, mx=mx)
        for (i in 1:nbpars)
            inits[i] <- dlist$inv.transforms[[i]](inits[i])
        res <- matrix(inits, ncol=1)
        dimnames(res) <- list(names(inits), "est")
        ret <- list(call=call, dlist=dlist, res=res, npars=0,
                    loglik=-minusloglik, AIC=2*minusloglik,
                    data = dat, datameans = colMeans(dat$X),
                    N=nrow(dat$Y), events=sum(dat$Y[,"status"]), trisk=sum(dat$Y[,"time"]))
    }
    else {
        optpars <- inits[setdiff(1:npars, fixedpars)]
        optim.args <- list(...)
        if (is.null(optim.args$method))
            optim.args$method <- "BFGS"
        gr <- if (!is.null(dlist$deriv)) Dminusloglik.flexsurv else NULL
        optim.args <- c(optim.args, list(par=optpars, fn=minusloglik.flexsurv, gr=gr,
                                         Y=Y, X=X, weights=weights, dlist=dlist,
                                         inits=inits, mx=mx, fixedpars=fixedpars,
                                         hessian=TRUE))
        opt <- do.call("optim", optim.args)
        est <- opt$par
        if (all(!is.na(opt$hessian)) && all(!is.nan(opt$hessian)) && all(is.finite(opt$hessian)) &&
            all(eigen(opt$hessian)$values > 0))
        {
            cov <- solve(opt$hessian); se <- sqrt(diag(cov))
            if (!is.numeric(cl) || length(cl)>1 || !(cl>0) || !(cl<1))
                stop("cl must be a number in (0,1)")
            lcl <- est - qnorm(1 - (1-cl)/2)*se
            ucl <- est + qnorm(1 - (1-cl)/2)*se
        }
        else {
            warning("Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite. ")
            cov <- lcl <- ucl <- NA
        }
        res <- cbind(est=inits, lcl=NA, ucl=NA)
        res[setdiff(1:npars, fixedpars),] <- cbind(est, lcl, ucl)
        colnames(res) <- c("est", paste(c("L","U"), round(cl*100), "%", sep=""))
        res.t <- res # results on transformed (log) scale
        for (i in 1:nbpars) # results on natural scale
            res[i,] <- dlist$inv.transforms[[i]](res[i,])
        ret <- list(call=match.call(), dlist=dlist, res=res, res.t=res.t, cov=cov,
                    coefficients=res.t[,"est"],
                    npars=length(est), fixedpars=fixedpars, optpars=setdiff(1:npars, fixedpars),
                    mx=mx, ncovs=ncovs, ncoveffs=ncoveffs, basepars=1:nbpars, covpars=(nbpars+1):npars,
                    loglik=-opt$value, AIC=2*opt$value + 2*length(est), cl=cl, opt=opt,
                    data = dat, datameans = colMeans(dat$X), terms=Terms,
                    N=nrow(dat$Y), events=sum(dat$Y[,"status"]), trisk=sum(dat$Y[,"time"]))
    }
    class(ret) <- "flexsurvreg"
    ret
}

### Compute CIs for survival, cumulative hazard and hazard at supplied
### times t and covariates X, using random sample of size B from the
### assumed MVN distribution of MLEs.

cisumm.flexsurvreg <- function(x, t, start, X, B=1000, cl=0.95) {
    if (any(is.na(x$res[,2])) || (B==0)) {
        ret <- array(NA, dim=c(length(t), 2, 3))
    }
    else {
        ret <- array(dim=c(B, length(t), 3))
        sim <- matrix(nrow=B, ncol=nrow(x$res))
        colnames(sim) <- rownames(x$res)
        sim[,x$optpars] <- rmvnorm(B, x$opt$par, x$cov)
        sim[,x$fixedpars] <- rep(x$res.t[x$fixedpars,"est"],each=B)
        for (i in seq(length=B)) {
            pcall <- list(q=t)
            for (j in x$dlist$pars)
                pcall[[j]] <- sim[i,j]
            beta <- if (x$ncoveffs==0) 0 else sim[i, x$covpars]
            for (j in seq(along=x$dlist$pars)){
                pcall[[x$dlist$pars[j]]] <- pcall[[x$dlist$pars[j]]] + X[x$mx[[j]]] %*% beta[x$mx[[j]]]
                pcall[[x$dlist$pars[j]]] <- x$dlist$inv.transforms[[j]](pcall[[x$dlist$pars[j]]])
            }
            probfn <- paste("p",x$dlist$name,sep="")
            surv <- 1 - do.call(probfn, pcall)
            dcall <- tcall <- pcall
            tcall$q <- start
            pobs <- 1 - do.call(probfn, tcall)
            densfn <- paste("d",x$dlist$name,sep="")
            names(dcall)[names(dcall)=="q"] <- "x"
            dens <- do.call(densfn, dcall)

            surv[t<start] <- 1; dens[t<start] <- 0
            ret[i,,] <- cbind(surv=surv/pobs, cumhaz=-log(surv/pobs), haz=pobs*dens/surv)
        }
        ret <- apply(ret, c(2,3), function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
        ret <- aperm(ret, c(2,1,3))
    }
    dimnames(ret)[[3]] <- c("surv","cumhaz","haz")
    ret
}

### Compute CIs for survival, cumulative hazard and hazard under
### spline-based models at supplied times t and covariates X, using
### random sample of size B from the assumed MVN distribution of MLEs.

cisumm.spline <- function(x, t, start, X, B=1000, cl=0.95) {
    if (any(is.na(x$res[,2])) || (B==0)) {
        ret <- array(NA, dim=c(length(t), 2, 3))
    }
    else {
        ret <- array(dim=c(B, length(t), 3))
        sim <- matrix(nrow=B, ncol=nrow(x$res))
        colnames(sim) <- rownames(x$res)
        sim[,x$optpars] <- rmvnorm(B, x$opt$par, x$cov)
        sim[,x$fixedpars] <- rep(x$res.t[x$fixedpars,"est"],each=B)
        for (i in seq(length=B)) {

            gamma <- sim[i, 1:(x$k + 2)]
            beta <- if (x$ncovs==0) 0 else sim[i, (x$k+3):(x$k + 2 + x$ncoveffs)]
            dens <- dsurvspline(t, gamma, beta, X, x$knots, x$scale)
            surv <- 1 - psurvspline(t, gamma, beta, X, x$knots, x$scale)
            pobs <- 1 - psurvspline(start, gamma, beta, X, x$knots, x$scale)

            surv[t<start] <- 1; dens[t<start] <- 0
            ret[i,,] <- cbind(surv=surv/pobs, cumhaz=-log(surv/pobs), haz=pobs*dens/surv)
        }
        ret <- apply(ret, c(2,3), function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
        ret <- aperm(ret, c(2,1,3))
    }
    dimnames(ret)[[3]] <- c("surv","cumhaz","haz")
    ret
}

print.flexsurvreg <- function(x, ...)
{
    covmeans <- colMeans(x$data$X)
    covs <- names(covmeans)
    covinds <- match(covs, rownames(x$res))
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (x$npars > 0) {
        res <- signif(x$res, 3)
        cat ("Estimates: \n")
        if (any(covinds)) {
            ecoefs <- matrix(NA, nrow=nrow(x$res), ncol=3)
            colnames(ecoefs) <- c("exp(est)", colnames(res)[2:3])
            means <- rep(NA,nrow(x$res))
            ecoefs[covinds,] <- signif(exp(x$res[covinds,1:3,drop=FALSE]), 3)
            means[covinds] <- signif(covmeans, 3)
            res <- cbind(means, res, ecoefs)
            colnames(res)[1] <- "data mean"
        }
        print(format(res), print.gap=2, quote=FALSE, na.print="")
    }
    cat("\nN = ", x$N, ",  Events: ", x$events,
        ",  Censored: ", x$N - x$events,
        "\nTotal time at risk: ", x$trisk,
        "\nLog-likelihood = ", x$loglik, ", df = ", x$npars,
        "\nAIC = ", x$AIC, "\n\n", sep="")
}

summary.flexsurvreg <- function(object, X=NULL, type="survival", t=NULL, start=NULL, B=1000, cl=0.95, ...)
{
    x <- object
    dat <- x$data
    isfac <- sapply(dat$Xraw,is.factor)
    ncovs <- x$ncovs
    type <- match.arg(type, c("survival","cumhaz","hazard"))
    if (ncovs > 0 && is.null(X)) {
        ## if any continuous covariates, calculate fitted survival for "average" covariate value
        if (!all(isfac))
            X <- matrix(colMeans(dat$X) ,nrow=1)
        ## else calculate for all different factor groupings
        else {
            X <- unique(dat$X)
            ## build names like "COVA=value1,COVB=value2"
            nam <- as.matrix(unique(dat$Xraw))
            for (i in 1:ncol(nam)) nam[,i] <- paste(colnames(nam)[i], nam[,i], sep="=")
            rownames(X) <- apply(nam, 1, paste, collapse=",")
        }
    }
    else if (is.null(X)) X <- as.matrix(0, nrow=1, ncol=max(ncol(dat$X),1))
    if (is.null(t)) {
        t <- sort(unique(dat$Y[,"stop"]))
        if (is.null(start))
            start <- dat$Y[order(dat$Y[!duplicated(dat$Y[,"stop"]),"stop"]),"start"]
    }
    else {
        if (is.null(start))
            start <- rep(0, length(t))
        else if (length(start) != length(t)) stop("length of \"t\" is ",length(t)," but length of \"start\" is ",length(start))
    }
    pcall <- list(q=t)
    if (!is.null(x$knots)) {
        gamma <- x$res[1:(x$k + 2),"est"]
        segamma <- x$res[1:(x$k + 2),"se"]
        beta <- if (ncovs==0) 0 else x$res[(x$k+3):(x$k + 2 + ncol(X)),"est"]
        sebeta <- if (ncovs==0) 0 else x$res[(x$k+3):(x$k + 2 + ncol(X)),"se"]
    }
    else beta <- if (ncovs==0) 0 else x$res[setdiff(rownames(x$res), x$dlist$pars),"est"]
    if (ncol(X) != length(beta)){
        isare <- if(length(beta)==1) "is" else "are"
        plural <- if(ncol(X)==1) "" else "s"
        pluralc <- if(length(beta)==1) "" else "s"
        stop("Supplied X has ", ncol(X), " column",plural," but there ",isare," ",
             length(beta), " covariate effect", pluralc)
    }
    dlist <- x$dlist
    ret <- vector(nrow(X), mode="list")
    names(ret) <- rownames(X)
    for (i in 1:nrow(X)) {
        if (is.null(x$knots)) {
            for (j in seq(along=dlist$pars)) {
                pcall[[dlist$pars[j]]] <- x$res[dlist$pars[j],"est"]
                mu <- dlist$transforms[[j]](pcall[[dlist$pars[j]]])
                mu <- mu + X[i,x$mx[[j]],drop=FALSE] %*% beta[x$mx[[j]]]
                pcall[[dlist$pars[j]]] <- dlist$inv.transforms[[j]](mu)
            }
            probfn <- paste("p",dlist$name,sep="")
            prob <- do.call(probfn, pcall)
            tcall <- pcall; tcall$q <- start
            pobs <- 1 - do.call(probfn, tcall) # =1 unless left-truncated
            res.ci <- cisumm.flexsurvreg(x, t, start, X[i,], B=B, cl=cl)
            prob[t<start] <- 0
            if (type=="survival") {
                y <- (1 - prob)/pobs
                ly <- res.ci[,1,"surv"]
                uy <-  res.ci[,2,"surv"]
            }
            else if (type=="cumhaz"){
                y <- -log((1 - prob)/pobs)
                ly <- res.ci[,1,"cumhaz"]
                uy <-  res.ci[,2,"cumhaz"]
            }
            else if (type=="hazard") {
                densfn <- paste("d",dlist$name,sep="")
                dcall <- pcall
                names(dcall)[names(dcall)=="q"] <- "x"
                dens <- do.call(densfn, dcall)
                dens[t<start] <- 0
                y <- pobs * dens / (1 - prob)
                ly <- res.ci[,1,"haz"]
                uy <-  res.ci[,2,"haz"]
            }
        }
        else {
            xd <- cbind(basis(x$knots, log(t)))
            nobs <- length(t)
            surv <- 1 - psurvspline(t, gamma, beta, X[i,], x$knots, x$scale)
            pobs <- 1 - psurvspline(start, gamma, beta, X[i,], x$knots, x$scale)
            res.ci <- cisumm.spline(x, t, start, X[i,], B=B, cl=cl)
            if (all(dat$Y[,"start"]==0)) {
                ## use analytic CIs if no left-truncation, else bootstrap-like CIs
                if (ncovs>0) xd <- cbind(xd, matrix(rep(X[i,],each=nobs),nrow=nobs))
                seeta <- numeric(nobs)
                for (j in 1:nobs) seeta[j] <- sqrt(xd[j,] %*% x$cov %*% xd[j,])
                lclsurv <- 1 - psurvspline(t, gamma, beta, X[i,], x$knots, x$scale, offset=qnorm(1 - (1-cl)/2)*seeta)
                uclsurv <- 1 - psurvspline(t, gamma, beta, X[i,], x$knots, x$scale, offset=-qnorm(1 - (1-cl)/2)*seeta)
                lclsurv[t==0] <- uclsurv[t==0] <- 1 # set by hand since basis() returns -Inf for log(0)
            }
            else { lclsurv <- res.ci[,1,"surv"]; uclsurv <- res.ci[,2,"surv"] }
            if (type=="survival") {y <- surv/pobs; ly <- lclsurv; uy <- uclsurv}
            else if (type=="cumhaz") {y <- -log(surv/pobs); ly <- -log(lclsurv); uy <- -log(uclsurv)}
            else if (type=="hazard") {
                dens <- dsurvspline(t, gamma, beta, X[i,], x$knots, x$scale)
                y <- dens/surv
                haz.ci <- cisumm.spline(x, t, start, X[i,], B=B, cl=cl)
                ly <- res.ci[,1,"haz"]
                uy <-  res.ci[,2,"haz"]
            }
        }
        ret[[i]] <- data.frame(time=t, est=y, lcl=ly, ucl=uy)
    }
    if (ncovs>0) ret$X <- X
    ret
}

plot.flexsurvreg <- function(x, X=NULL, type="survival", t=NULL, start=NULL,
                             est=TRUE, ci=NULL, B=1000, cl=0.95,
                             col.obs="black", lty.obs=1, lwd.obs=1,
                             col="red",lty=1,lwd=2,
                             col.ci=NULL,lty.ci=2,lwd.ci=1,
                             add=FALSE,...)
{
    type <- match.arg(type, c("survival","cumhaz","hazard"))
    ## don't calculate or plot CIs by default if all covs are categorical -> multiple curves
    if (is.null(ci))
        ci <- ((x$ncovs == 0) || (!(sapply(x$data$Xraw,is.factor))))
    if (!ci) B <- 0
    summ <- summary.flexsurvreg(x, X=X, type=type, t=t, B=B, cl=cl)
    t <- summ[[1]]$time
    X <- if (is.null(summ$X)) as.matrix(0, nrow=1, ncol=max(x$ncoveffs,1)) else summ$X
    if (is.null(col.ci)) col.ci <- col
    if (is.null(lwd.ci)) lwd.ci <- lwd
    dat <- x$data
    ncovs <- x$ncovs
    isfac <- sapply(dat$Xraw,is.factor)
    if (!add) {
        form <- "Surv(dat$Y[,\"start\"],dat$Y[,\"stop\"],dat$Y[,\"status\"]) ~ "
        form <- paste(form, if (ncovs > 0 && all(isfac)) paste("dat$X[,",1:x$ncoveffs,"]", collapse=" + ") else 1)
        form <- as.formula(form)
        ## If any continuous covariates, it is hard to define subgroups
        ## so just plot the population survival
        if (type=="survival") {
            plot(survfit(form, data=as.data.frame(dat$X)), col=col.obs, lty=lty.obs, lwd=lwd.obs, ...)
        }
        else if (type=="cumhaz") {
            plot(survfit(form, data=as.data.frame(dat$X)), fun="cumhaz", col=col.obs, lty=lty.obs, lwd=lwd.obs, ...)
        }
        else if (type=="hazard") {
            if (!all(dat$Y[,"start"]==0)) warning("Left-truncated data not supported by muhaz: ignoring truncation point when plotting observed hazard")
            if (!all(isfac))
                plot(muhaz(dat$Y[,"stop"], dat$Y[,"status"], ...),
                     col=col.obs, lty=lty.obs, lwd=lwd.obs)
            else {
                ## plot hazard for all groups defined by unique combinations of covariates
                group <- if(ncovs>0) do.call("interaction", as.data.frame(dat$X)) else factor(rep(1,nrow(dat$Y)))
                haz <- list()
                for (i in 1:nrow(X)) {
                    subset <- (group == unique(group)[i])
                    haz[[i]] <- muhaz(dat$Y[,"time"], dat$Y[,"status"], subset=subset, ...)
                }
                plot(haz[[1]], col=col.obs, lty=lty.obs, lwd=lwd.obs,
                     ylim=range(sapply(haz, function(x)range(x$haz.est))), ...)
                if (nrow(X)>1) {
                    for (i in 1:nrow(X)) {
                        lines(haz[[i]], col=col.obs, lty=lty.obs, lwd=lwd.obs)
                    }
                }
            }
        }
    }
    for (i in 1:nrow(X)) {
        if (est) lines(summ[[i]]$t, summ[[i]]$est, col=col, lty=lty, lwd=lwd)
        if (ci) {
            lines(summ[[i]]$t, summ[[i]]$lcl, col=col.ci, lty=lty.ci, lwd=lwd.ci)
            lines(summ[[i]]$t, summ[[i]]$ucl, col=col.ci, lty=lty.ci, lwd=lwd.ci)
        }
    }
}

lines.flexsurvreg <- function(x, X=NULL, type="survival", t=NULL,
                              est=TRUE, ci=NULL, B=1000, cl=0.95,
                              col="red",lty=1,lwd=2,
                              col.ci=NULL,lty.ci=2,lwd.ci=1, ...)
{
    plot.flexsurvreg(x, X, type=type, t=t, est=est, ci=ci, B=B, cl=cl,
                     col=col, lty=lty, lwd=lwd, col.ci=col.ci,lty.ci=lty.ci,lwd.ci=lwd.ci, add=TRUE, ...)
}

vcov.flexsurvreg <- function (object, ...)
{
    object$cov
}
