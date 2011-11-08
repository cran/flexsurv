
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
        ret <- 1 / shape * log1p(-log1p(-p) * shape / rate)
    }
    ret
}

rgompertz <- function(n, shape = 1, rate = 1){
    if (!check.gompertz(shape=shape, rate=rate)) return(rep(NaN, n))
    check.gompertz(shape=shape, rate=rate)
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
    if (rate<=0) {warning("Non-positive rate parameter"); ret <- FALSE}
    ret
}
