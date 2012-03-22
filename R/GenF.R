## Generalized F distribution (Prentice 1975 parameterisation) 
## For P=0 this is equivalent to the generalized (log-)gamma (Prentice 1974) 
## P=Q=0, lognormal
## P=0, Q=1, Weibull 
## Equation 2 in C.Cox (2008) is wrong, delta*beta*m1 not beta*m1 in first exponent in numerator

dgenf <- function(x, mu=0, sigma=1, Q, P, log=FALSE) { 
    if (!check.genf(mu=mu, sigma=sigma, Q=Q, P=P)) return(rep(NaN, length(x)))
    ret <- numeric(length(x))
    ret[x<=0] <- 0
    xx <- x[x>0]
    if (P==0) {
        logdens <- dgengamma(xx, mu, sigma, Q, log=TRUE)
    }
    else {
        tmp <- Q^2 + 2*P
        delta <- sqrt(tmp)
        s1 <- 2 / (tmp + Q*delta)
        s2 <- 2 / (tmp - Q*delta)
        expw <- xx^(delta/sigma)*exp(-mu*delta/sigma)
        logdens <- log(delta) + s1/sigma*delta*(log(xx) - mu) + s1*(log(s1) - log(s2)) -
            log(sigma*xx) - (s1+s2)*log(1 + s1*expw/s2) - lbeta(s1, s2)
    }
    ret[x>0] <- if (log) logdens else exp(logdens)
    ret
}

pgenf <- function(q, mu=0, sigma=1, Q, P, lower.tail = TRUE, log.p = FALSE)
{
    if (!check.genf(mu=mu, sigma=sigma, Q=Q, P=P)) return(rep(NaN, length(q)))
    q[q<0] <- 0
    if (P==0) {
        ret <- pgengamma(q, mu, sigma, Q, lower.tail=TRUE, log.p=FALSE)
    }
    else {
        tmp <- Q^2 + 2*P
        delta <- sqrt(tmp)
        s1 <- 2 / (tmp + Q*delta)
        s2 <- 2 / (tmp - Q*delta)
        expw <- q^(delta/sigma)*exp(-mu*delta/sigma)    
        ret <- 1 - pbeta(s2/(s2 + s1*expw), s2, s1)
    }
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

Hgenf <- function(x, mu=0, sigma=1, Q, P)
{
    -log(pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail=FALSE))
}

hgenf <- function(x, mu=0, sigma=1, Q, P)
{
    dgenf(x=x, mu=mu, sigma=sigma, Q=Q, P=P) / 
        pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail=FALSE)
}

qgenf <- function(p, mu=0, sigma=1, Q, P, lower.tail = TRUE, log.p = FALSE)
{
    if (!check.genf(mu=mu, sigma=sigma, Q=Q, P=P)) return(rep(NaN, length(p)))
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    if (P==0) {
        ret <- qgengamma(p, mu, sigma, Q)
    }
    else { 
        tmp <- Q^2 + 2*P
        delta <- sqrt(tmp)
        s1 <- 2 / (tmp + Q*delta)
        s2 <- 2 / (tmp - Q*delta)
        w <- log(qf(p, 2*s1, 2*s2))
        ret <- exp(w*sigma/delta + mu)
    }
    ret
}

rgenf <- function(n, mu=0, sigma=1, Q, P)
{
    if (!check.genf(mu=mu, sigma=sigma, Q=Q, P=P)) return(rep(NaN, n))
    if (P==0) {
        ret <- rgengamma(n, mu, sigma, Q)
    }
    else { 
        tmp <- Q^2 + 2*P
        delta <- sqrt(tmp)
        s1 <- 2 / (tmp + Q*delta)
        s2 <- 2 / (tmp - Q*delta)
        w <- log(rf(n, 2*s1, 2*s2))
        ret <- exp(w*sigma/delta + mu)
    }
    ret
}

check.genf <- function(mu, sigma, Q, P){
    ret <- TRUE
    if (missing(Q)) stop("shape parameter \"Q\" not given")
    if (missing(P)) stop("shape parameter \"P\" not given")
    if (any(sigma <= 0)) {warning("Non-positive scale parameter \"sigma\""); ret <- FALSE}
    if (any(P < 0)) {warning("Negative shape parameter \"P\""); ret <- FALSE}
    ret
}

dgenf.orig <- function(x, mu=0, sigma=1, s1, s2, log=FALSE) { 
    if (!check.genf.orig(mu=mu, sigma=sigma, s1=s1, s2=s2)) return(rep(NaN, length(x)))
    ret <- numeric(length(x))
    ret[x<=0] <- 0
    xx <- x[x>0]
    w <- (log(xx) - mu)/sigma
    expw <- xx^(1/sigma)*exp(-mu/sigma)
    logdens <- -log(sigma*xx) + s1*(log(s1) + w - log(s2)) - (s1+s2)*log(1 + s1*expw/s2) - lbeta(s1, s2)
    ret[x>0] <- if (log) logdens else exp(logdens)
    ret
}

pgenf.orig <- function(q, mu=0, sigma=1, s1, s2, lower.tail = TRUE, log.p = FALSE)
{
    if (!check.genf.orig(mu=mu, sigma=sigma, s1=s1, s2=s2)) return(rep(NaN, length(q)))
    q[q<0] <- 0
    w <- (log(q) - mu)/sigma    
    ret <- 1 - pbeta(s2/(s2 + s1*exp(w)), s2, s1)
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

Hgenf.orig <- function(x, mu=0, sigma=1, s1, s2)
{
    -log(pgenf.orig(q=x, mu=mu, sigma=sigma, s1=s1, s2=s2, lower.tail=FALSE))
}

hgenf.orig <- function(x, mu=0, sigma=1, s1, s2)
{
    dgenf.orig(x=x, mu=mu, sigma=sigma, s1=s1, s2=s2) / 
        pgenf.orig(q=x, mu=mu, sigma=sigma, s1=s1, s2=s2, lower.tail=FALSE)
}

qgenf.orig <- function(p, mu=0, sigma=1, s1, s2, lower.tail = TRUE, log.p = FALSE)
{
    if (!check.genf.orig(mu=mu, sigma=sigma, s1=s1, s2=s2)) return(rep(NaN, length(p)))
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    w <- log(qf(p, 2*s1, 2*s2))
    exp(w*sigma + mu)
}

rgenf.orig <- function(n, mu=0, sigma=1, s1, s2)
{
    if (!check.genf.orig(mu=mu, sigma=sigma, s1=s1, s2=s2)) return(rep(NaN, n))
    w <- log(rf(n, 2*s1, 2*s2))
    exp(w*sigma + mu)
}

check.genf.orig <- function(mu, sigma, s1, s2){
    ret <- TRUE
    if (missing(s1)) stop("shape parameter \"s1\" not given")
    if (missing(s2)) stop("shape parameter \"s2\" not given")
    if (any(sigma <= 0)) {warning("Non-positive scale parameter \"sigma\""); ret <- FALSE}
    if (any(s1 <= 0)) {warning("Non-positive shape parameter \"s1\""); ret <- FALSE}
    if (any(s2 <= 0)) {warning("Non-positive shape parameter \"s2\""); ret <- FALSE}
    ret
}

## Thanks to Skif Pankov
mean.genf.orig <- function(mu, sigma, s1, s2){
    exp(mu) * (s2/s1)^sigma * gamma(s1 + sigma)*gamma(s2 - sigma) / (gamma(s1)*gamma(s2))
}

