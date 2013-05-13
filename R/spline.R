flexsurvspline <- function(formula, data, k=0, knots=NULL, scale="hazard", inits=NULL, fixedpars=NULL, cl=0.95, ...)
{
    call <- match.call()
    indx <- match(c("formula", "data"), names(call), nomatch = 0)
    if (indx[1] == 0)
        stop("A \"formula\" argument is required")
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    m <- eval(temp, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    if (!(attr(Y, "type")  %in% c("right","counting")))
        stop("Survival object type \"", attr(Y, "type"), "\"", " not supported")
    if (attr(Y, "type") == "counting")
        Y <- cbind(Y, time=Y[,"stop"] - Y[,"start"]) # converts Y from Surv object to numeric matrix
    else Y <- cbind(Y, start=0, stop=Y[,"time"])
    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m)
    dat <- list(Y=Y, X=X[,-1,drop=FALSE], Xraw=m[,-1,drop=FALSE])
    X <- dat$X
    ## knots=m=0, df=1, == no knots = weibull
    ## choose knot locations by quantiles
    ## inits : m+2 gammas, beta
    if (is.null(knots)) {
        is.wholenumber <-
            function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
        if (is.null(k)) stop("either \"knots\" or \"k\" must be specified")
        if (!is.numeric(k)) stop("k must be numeric")
        if (!is.wholenumber(k) || (k<0)) stop("number of knots \"k\" must be a non-negative integer")
        knots <- quantile(log(Y[,"stop"])[Y[,"status"]==1], seq(0, 1, length=k+2))
    }
    else {
        if (!is.numeric(knots)) stop("\"knots\" must be a numeric vector")
        minlogtime <- min(log(Y[,"stop"]))
        if (any(knots <= minlogtime)) {
            badknots <- knots[knots < min(log(Y[,"stop"]))]
            plural <- if (length(badknots) > 1) "s" else ""
            stop("knot", plural, " ", paste(badknots,collapse=", "), " less than or equal to minimum log time ", minlogtime)
        }
        maxlogtime <- max(log(Y[,"stop"]))
        if (any(knots >= maxlogtime)) {
            badknots <- knots[knots > maxlogtime]
            plural <- if (length(badknots) > 1) "s" else ""
            stop("knot", plural, " ", paste(badknots,collapse=", "), " greater than or equal to maximum log time ", maxlogtime)
        }
        k <- length(knots)
        knots <- c(min(log(Y[,"stop"])[Y[,"status"]==1]), knots, max(log(Y[,"stop"])[Y[,"status"]==1]))
    }
    match.arg(scale, c("hazard","odds","normal"))
    ncovs <- ncol(dat$Xraw)
    ncoveffs <- ncol(X)
    npars <- k + 2 + ncoveffs
    if (is.null(inits)) {
        inits <- flexsurv.splineinits(Y, X, data, knots, scale)
    }
    else {
        if (!is.numeric(inits) || length(inits) != k + 2 + ncoveffs)
            stop("inits must be a numeric vector of length ", k + 2 + ncoveffs, ": 2 + ",
                 k, " knots + ", ncoveffs, " covariate effects")
    }
    cnames <- if(ncoveffs==0) NULL else colnames(X)
    names(inits) <- c(paste("gamma",0:(k+1),sep=""), cnames)
    if (!is.null(fixedpars) && !is.logical(fixedpars) &&
        (!is.numeric(fixedpars) || any(!(fixedpars %in% 1:npars)))){
        dots <- if(npars>2) "...," else ""
        stop("fixedpars must be TRUE/FALSE or a vector of indices in 1,",dots,npars)
    }
    if ((is.logical(fixedpars) && fixedpars==TRUE) ||
        (is.numeric(fixedpars) && all(fixedpars == 1:npars))) {
        minusloglik <- minusloglik.stpm(optpars=inits, knots=knots, Y=Y, X=X, inits=inits,
                                   fixedpars=NULL, scale=scale)
        res <- cbind(est=inits,lcl=NA,ucl=NA)
        ret <- list(call=match.call(), k=k, knots=knots, scale=scale, res=res, npars=0,
                    loglik=-minusloglik, AIC=2*minusloglik,
                    data = dat, datameans = colMeans(dat$X),
                    N=nrow(dat$Y), events=sum(dat$Y[,"status"]), trisk=sum(dat$Y[,"time"]))
    }
    else {
        optpars <- inits[setdiff(1:npars, fixedpars)]
        optim.args <- list(...)
        if (is.null(optim.args$method))
            optim.args$method <- "BFGS"
        gr <- if (scale=="normal") NULL else Dminusloglik.stpm
        optim.args <- c(optim.args, list(par=optpars, fn=minusloglik.stpm, gr=gr,
                                         knots=knots, Y=Y, X=X,
                                         inits=inits, fixedpars=fixedpars,
                                         scale=scale, hessian=TRUE))
        opt <- do.call("optim", optim.args)
        est <- opt$par
        cov <- solve(opt$hessian); se <- sqrt(diag(cov))
        if (!is.numeric(cl) || length(cl)>1 || !(cl>0) || !(cl<1))
            stop("cl must be a number in (0,1)")
        lcl <- est - qnorm(1 - (1-cl)/2)*se
        ucl <- est + qnorm(1 - (1-cl)/2)*se
        res <- cbind(est=inits, lcl=NA, ucl=NA, se=NA)
        res[setdiff(1:npars, fixedpars),] <- cbind(est, lcl, ucl, se)
        colnames(res) <- c("est", paste(c("L","U"), round(cl*100), "%", sep=""), "se")
        ret <- list(call=match.call(), k=k, knots=knots, scale=scale, res=res, cov=cov,
                    npars=length(est), fixedpars=fixedpars, optpars=setdiff(1:npars, fixedpars),
                    ncovs=ncovs, ncoveffs=ncoveffs,
                    loglik=-opt$value, AIC=2*opt$value + 2*length(est), cl=cl, opt=opt,
                    data = dat, datameans = colMeans(dat$X),
                    N=nrow(dat$Y), events=sum(dat$Y[,"status"]), trisk=sum(dat$Y[,"time"]))
    }
    class(ret) <- "flexsurvreg"
    ret
}

flexsurv.splineinits <- function(Y, X, data, knots, scale)
{
    X <- X[Y[,"status"]==1,,drop=FALSE]
    Y <- Y[Y[,"status"]==1,,drop=FALSE]
    ## using coxph on original formula followed by survfit.coxph fails
    ## due to scoping
    form <- paste("Surv(Y[,\"time\"], Y[,\"status\"]) ~ ")
    if (ncol(X)>0)
        form <- paste(form, paste(paste("X[,",1:ncol(X),"]",sep=""), collapse=" + "))
    else form <- paste(form, "1")
    cox <- coxph(as.formula(form))
    surv <- survfit(cox, data=cbind(Y, X))
    surv <- surv$surv[match(Y[,"time"], surv$time)]
    if (scale=="hazard")
        logH <- log(-log(surv))
    else if (scale=="odds")
        logH <- log((1 - surv)/surv)
    else if (scale=="normal")
        logH <- qnorm(1 - surv)
    b <- basis(knots, log(Y[,"time"]))
    form <- paste("logH ~ ",
                  paste(paste("b[,",2:ncol(b),"]",sep=""), collapse=" + "))
    if (ncol(X)>0)
        form <- paste(form, "+", paste(paste("X[,",1:ncol(X),"]",sep=""), collapse=" + "))
    inits <- coef(lm(as.formula(form)))
    inits
}

minusloglik.stpm <- function(optpars, knots, Y, X=0, inits, fixedpars=NULL, scale="hazard"){
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
    dens <- dsurvspline(Y[dead,"stop"], gamma, beta, X[dead,,drop=FALSE], knots, scale)
    surv <- 1 - psurvspline(Y[!dead,"stop"], gamma, beta, X[!dead,,drop=FALSE], knots, scale)
    pobs <- 1 - psurvspline(Y[,"start"], gamma, beta, X[,,drop=FALSE], knots, scale) # = 1 unless left-truncated
    ## workaround to avoid warnings, TODO think about implicit parameter constraints instead
    if (any(dens<=0) || any(surv<=0)) return(Inf)
    - ( sum(log(dens)) + sum(log(surv)) - sum(log(pobs)))
}

psurvspline <- function(q, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0){
    if (length(gamma) != length(knots)) stop("length of gamma should equal number of knots")
    match.arg(scale, c("hazard","odds","normal"))
    eta <- basis(knots, log(q)) %*% gamma + as.numeric(X %*% beta) + offset
    surv <- if (scale=="hazard")  exp(-exp(eta)) else if (scale=="odds") 1 / (1 + exp(eta)) else if (scale=="normal") pnorm(-eta)
    as.numeric(1 - surv)
}

dsurvspline <- function(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0){
    if (length(gamma) != length(knots)) stop("length of gamma should equal number of knots")
    match.arg(scale, c("hazard","odds","normal"))
    eta <- basis(knots, log(x)) %*% gamma + as.numeric(X %*% beta) + offset
    eeta <- if (scale=="hazard") exp(eta - exp(eta)) else if (scale=="odds")  exp(eta) / (1 + exp(eta))^2 else if (scale=="normal") dnorm(eta)
    dens <- 1 / x * dbasis(knots, log(x)) %*% gamma * eeta
    as.numeric(dens)
}

basis <- function(knots, x) {
    nk <- length(knots)
    lam <- (knots[nk] - knots)/(knots[nk] - knots[1])
    b <- matrix(nrow=length(x), ncol=nk)
    b[,1] <- 1
    b[,2] <- x
    if (nk>2) {
        for (j in 1:(nk-2)) {
            b[,j+2] <- pmax(x - knots[j+1], 0)^3 - lam[j+1]*pmax(x - knots[1], 0)^3 -
                (1 - lam[j+1])*pmax(x - knots[nk], 0)^3
        }
    }
    b
}

dbasis <- function(knots, x) {
    nk <- length(knots)
    lam <- (knots[nk] - knots)/(knots[nk] - knots[1])
    b <- matrix(nrow=length(x), ncol=nk)
    b[,1] <- 0
    b[,2] <- 1
    if (nk>2) {
        for (j in 3:nk) {
            b[,j] <- 3*pmax(x - knots[j-1], 0)^2 - 3*lam[j-1]*pmax(x - knots[1], 0)^2 -
                3*(1 - lam[j-1])*pmax(x - knots[nk], 0)^2
        }
    }
    b
}
