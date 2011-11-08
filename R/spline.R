flexsurvspline <- function(formula, data, k=0, knots=NULL, scale="hazard", inits=NULL, fixedpars=NULL, cl=0.95,...)
{
    call <- match.call()
    dat <- form.survdata(call, formula, data)
    Y <- dat$Y; X <- dat$X 
    ## knots=m=0, df=1, == no knots = weibull 
    ## choose knot locations by quantiles
    ## inits : m+2 gammas, beta 
    if (is.null(knots)) {
        is.wholenumber <-
            function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
        if (!is.wholenumber(k) || (k<0)) stop("number of knots \"k\" must be a non-negative integer")
        knots <- quantile(log(Y[,"time"])[Y[,"status"]==1], seq(0, 1, length=k+2))
    }
    else {
        if (!is.numeric(knots)) stop("\"knots\" must be a numeric vector")
        minlogtime <- min(log(Y[,"time"]))
        if (any(knots <= minlogtime)) {
            badknots <- knots[knots < min(log(Y[,"time"]))]
            plural <- if (length(badknots) > 1) "s" else ""
            stop("knot", plural, " ", paste(badknots,collapse=", "), " less than or equal to minimum log time ", minlogtime)
        }
        maxlogtime <- max(log(Y[,"time"]))
        if (any(knots >= maxlogtime)) {
            badknots <- knots[knots > maxlogtime]
            plural <- if (length(badknots) > 1) "s" else ""
            stop("knot", plural, " ", paste(badknots,collapse=", "), " greater than or equal to maximum log time ", maxlogtime)
        }
        k <- length(knots)
        knots <- c(minlogtime, knots, maxlogtime)
    } 
    match.arg(scale, c("hazard","odds","normal"))
    ncovs <- ncol(X)
    npars <- k + 2 + ncovs
    if (is.null(inits)) {
        inits <- flexsurv.splineinits(Y, X, data, knots, scale)
    }
    else { 
        if (!is.numeric(inits) || length(inits) != k + 2 + ncovs)
            stop("inits must be a numeric vector of length ", k + 2 + ncovs, ": 2 + ",
                 k, " knots + ", ncovs, " covariates")
    }
    cnames <- if(ncovs==0) NULL else colnames(X)
    names(inits) <- c(paste("gamma",0:(k+1),sep=""), cnames)
    if (!is.null(fixedpars) && !is.logical(fixedpars) &&
        (!is.numeric(fixedpars) || any(!(fixedpars %in% 1:npars)))){
        dots <- if(npars>2) "...," else ""
        stop("fixedpars must be TRUE/FALSE or a vector of indices in 1,",dots,npars)
    }
    if ((is.logical(fixedpars) && fixedpars==TRUE) ||
        (is.numeric(fixedpars) && all(fixedpars == 1:npars))) {
        minusloglik <- minusloglik.stpm(optpars=inits, knots=knots, t=Y[,"time"], dead=Y[,"status"], X=X, inits=inits,
                                   fixedpars=NULL, scale=scale)
        res <- cbind(est=inits,lcl=NA,ucl=NA)
        ret <- list(call=match.call(), k=k, knots=knots, scale=scale, res=res,
                    loglik=-minusloglik, AIC=2*minusloglik + 2*npars)
    }
    else { 
        optpars <- inits[setdiff(1:npars, fixedpars)]
        opt <- optim(optpars, minusloglik.stpm, knots=knots,
                     t=Y[,"time"], dead=Y[,"status"], X=X,
                     inits=inits, fixedpars=fixedpars, scale=scale, hessian=TRUE, ...)
        est <- opt$par
        cov <- solve(opt$hessian); se <- sqrt(diag(cov))
        if (!is.numeric(cl) || length(cl)>1 || !(cl>0) || !(cl<1))
            stop("cl must be a number in (0,1)")
        lcl <- est - qnorm(1 - (1-cl)/2)*se
        ucl <- est + qnorm(1 - (1-cl)/2)*se
        res <- cbind(est=inits, lcl=NA, ucl=NA)
        res[setdiff(1:npars, fixedpars),] <- cbind(est, lcl, ucl)
        colnames(res) <- c("est", paste(c("L","U"), round(cl*100), "%", sep=""))
        ret <- list(call=match.call(), k=k, knots=knots, scale=scale, res=res, 
                    loglik=-opt$value, AIC=2*opt$value + 2*npars)
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

minusloglik.stpm <- function(optpars, knots, t, dead, X=0, inits, fixedpars=NULL, scale="hazard"){
    pars <- inits
    npars <- length(pars)
    pars[setdiff(1:npars, fixedpars)] <- optpars
    x <- log(t)
    nk <- length(knots)
    gamma <- pars[1:nk] # always at least two of these, intercept + weibull par 
    eta <- fs.spline(gamma, x, knots)
    if (npars > nk) { # any covariates? 
        beta <- pars[(nk+1):npars]
        eta <- eta + X %*% beta
    }
    dead <- as.logical(dead)
    if (scale=="hazard") { 
        dens <- 1 / t[dead] * fs.dspline(gamma, x[dead], knots) * exp(eta[dead] - exp(eta[dead]))
        surv <- exp(-exp(eta[!dead]))
    }
    else if (scale=="odds") {
        dens <- 1 / t[dead] * fs.dspline(gamma, x[dead], knots) * exp(eta[dead]) / (1 + exp(eta[dead]))^2
        surv <- 1 / (1 + exp(eta[!dead]))
    }
    else if (scale=="normal") {
        dens <- 1 / t[dead] * fs.dspline(gamma, x[dead], knots) * dnorm(eta[dead])
        surv <- pnorm(-eta[!dead])
    }
    - ( sum(log(dens)) + sum(log(surv)) )
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

fs.spline <- function(gamma, x, knots){
    b <- basis(knots, x)
    spline <- b %*% gamma
    spline
} 

fs.dspline <- function(gamma, x, knots){
    nk <- length(knots)
    lam <- (knots[nk] - knots)/(knots[nk] - knots[1])
    basis <- matrix(nrow=length(x), ncol=nk-1)
    basis[,1] <- 1
    if (nk>2) { 
        for (j in 2:(nk-1)) {
            basis[,j] <- 3*pmax(x - knots[j], 0)^2 - 3*lam[j]*pmax(x - knots[1], 0)^2 -
                3*(1 - lam[j])*pmax(x - knots[nk], 0)^2
        }
    }    
    dspline <- basis %*% gamma[-1]
    dspline    
}
