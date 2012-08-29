
## consistent with paper and Stata 
## shape pos: inc haz, shape neg: dec haz just like weibull, shape zero=exponential   
## 
## in eha: shape=lam, gamma=1/scale
## log(shape) + x/scale - shape * scale * (exp(x/scale) - 1))
## shape/scale labelled wrong way round.

dgompertz <- function(x, shape, rate=1, log=FALSE) {
    if (!check.gompertz(shape=shape, rate=rate)) return(rep(NaN, length(x)))
    ret <- numeric(length(x))
    ret[x<0] <- 0
    t <- x[x>=0]
    if (shape==0) {
        logdens <- dexp(t, rate=rate, log=TRUE)
    }
    else { 
        logdens <- log(rate) + shape*t - rate/shape*(exp(shape*t) - 1)
    }    
    ret[x>=0] <- if (log) logdens else exp(logdens)
    ret
}

pgompertz <- function(q, shape, rate=1, lower.tail = TRUE, log.p = FALSE) {
    if (!check.gompertz(shape=shape, rate=rate)) return(rep(NaN, length(q)))
    q[q<0] <- 0
    if (shape==0) {
        ret <- pexp(q, rate=rate)
    }
    else { 
        ret <- 1 - exp(-rate/shape*(exp(shape*q) - 1))
    }
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

qgompertz <- function(p, shape, rate=1, lower.tail = TRUE, log.p = FALSE) {
    if (!check.gompertz(shape=shape, rate=rate)) return(rep(NaN, length(p)))
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    if (shape==0) {
        ret <- qexp(p, rate=rate)
    }
    else {
        asymp <- 1 - exp(rate/shape)
        immortal <- shape < 0 & p > asymp
        ret <- numeric(length(p))
        ret[immortal] <- Inf

        ## replicate arguments to length of longest vector (ugh)
        if (length(shape) < length(p))
            shape <- rep(shape, length=length(p))
        if (length(rate) < length(p))
            shape <- rep(shape, length=length(p))
        if (length(p) < length(rate))
            p <- rep(p, length=length(rate))
        if (length(p) < length(shape))
            p <- rep(p, length=length(shape))
        if (length(shape) < length(rate))
            shape <- rep(shape, length=length(rate))
        if (length(rate) < length(shape))
            rate <- rep(rate, length=length(shape))
        
        ret[!immortal] <- 1 / shape[!immortal] *
            log1p(-log1p(-p[!immortal]) * shape[!immortal] / rate[!immortal])
    }
    ret
}

rgompertz <- function(n, shape = 1, rate = 1){
    if (!check.gompertz(shape=shape, rate=rate)) return(rep(NaN, n))
    qgompertz(p=runif(n), shape=shape, rate=rate)
}

hgompertz <- function(x, shape, rate=1){
    rate*exp(shape*x)
}

Hgompertz <- function(x, shape, rate=1){
    if (shape==0) 
        rate*x
    else 
        rate/shape*expm1(shape*x)
}

check.gompertz <- function(shape, rate=1){
    ret <- TRUE
    if (missing(shape)) stop("shape parameter not given")
    if (any(rate<=0)) {warning("Non-positive rate parameter"); ret <- FALSE}
    ret
}
