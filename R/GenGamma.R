## Log-gamma or generalized gamma distribution (parameterisation as in Farewell and Prentice, Technometrics 1977)

dgengamma <- function(x, mu=0, sigma=1, Q, log=FALSE) {
    if (!check.gengamma(mu=mu, sigma=sigma, Q=Q)) return(rep(NaN, length(x)))
    ret <- numeric(length(x))
    ret[x<=0] <- 0
    xx <- x[x>0]
    if (Q != 0) { 
        y <- log(xx)
        w <- (y - mu)/sigma
        logdens <- -log(sigma*xx) + log(abs(Q)) + (Q^-2)*log(Q^-2) + Q^-2*(Q*w - exp(Q*w)) - lgamma(Q^-2)
    }
    else logdens <- dlnorm(xx, mu, sigma, log=TRUE)
    ret[x>0] <- if (log) logdens else exp(logdens)
    ret
}

pgengamma <- function(q, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE) {
    if (!check.gengamma(mu=mu, sigma=sigma, Q=Q)) return(rep(NaN, length(q)))
    q[q<0] <- 0
    if (Q != 0) { 
        y <- log(q)
        w <- (y - mu)/sigma
        expnu <- exp(Q*w)*Q^-2
        ret <- if (Q > 0) pgamma(expnu, Q^-2) else 1 - pgamma(expnu, Q^-2)
    }
    else ret <- plnorm(q, mu, sigma)
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

Hgengamma <- function(x, mu=0, sigma=1, Q)
{
    -log(pgengamma(q=x, mu=mu, sigma=sigma, Q=Q, lower.tail=FALSE))
}

hgengamma <- function(x, mu=0, sigma=1, Q)
{
    dgengamma(x=x, mu=mu, sigma=sigma, Q=Q) / 
        pgengamma(q=x, mu=mu, sigma=sigma, Q=Q, lower.tail=FALSE)
}

qgengamma <- function(p, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE)
{
    if (!check.gengamma(mu=mu, sigma=sigma, Q=Q)) return(rep(NaN, length(p)))
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    if (Q==0)
        qlnorm(p, mu, 1/sigma^2)
    else {
        if (Q < 0) p <- 1 - p
        w <- log(Q^2*qgamma(p, 1/Q^2, 1)) / Q
        exp(mu + sigma*w)
    }
}

rgengamma <- function(n, mu=0, sigma=1, Q) {
    if (!check.gengamma(mu=mu, sigma=sigma, Q=Q)) return(rep(NaN, n))
    check.gengamma(mu=mu, sigma=sigma, Q=Q)
    if (Q==0)
        rlnorm(n, mu, 1/sigma^2)
    else {
        w <- log(Q^2*rgamma(n, 1/Q^2, 1)) / Q
        exp(mu + sigma*w)
    }
}

check.gengamma <- function(mu, sigma, Q){
    ret <- TRUE
    if (missing(Q)) stop("shape parameter \"Q\" not given")
    if (any(sigma <= 0)) {warning("Non-positive scale parameter \"sigma\""); ret <- FALSE}
    ret
}

dgengamma.orig <- function(x, shape, scale=1, k, log=FALSE){
    if (!check.gengamma.orig(shape=shape, scale=scale, k=k)) return(rep(NaN, length(x)))
    ret <- numeric(length(x))
    ret[x<=0] <- 0
    t <- x[x>0]
    logdens <- log(shape) - lgamma(k) + (shape*k - 1)*log(t) - shape*k*log(scale) - (t/scale)^shape 
    ret[x>0] <- if (log) logdens else exp(logdens)
    ret    
}

pgengamma.orig <- function(q, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE) {
    if (!check.gengamma.orig(shape=shape, scale=scale, k=k)) return(rep(NaN, length(q)))
    q[q<0] <- 0
    y <- log(q)
    w <- (y - log(scale))*shape
    ret <- pgamma(exp(w), shape=k) 
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

Hgengamma.orig <- function(x, shape, scale=1, k)
{
    -log(pgengamma.orig(q=x, shape=shape, scale=scale, k=k, lower.tail=FALSE))
}

hgengamma.orig <- function(x, shape, scale=1, k)
{
    dgengamma.orig(x=x, shape=shape, scale=scale, k=k) / 
        pgengamma.orig(q=x, shape=shape, scale=scale, k=k, lower.tail=FALSE)
}

qgengamma.orig <- function(p, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE)
{
    if (!check.gengamma.orig(shape=shape, scale=scale, k=k)) return(rep(NaN, length(p)))
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    w <- log(qgamma(p, shape=k))
    y <- w / shape  + log(scale)
    exp(y) 
}

rgengamma.orig <- function(n, shape, scale=1, k) {
    if (!check.gengamma.orig(shape=shape, scale=scale, k=k)) return(rep(NaN, n))
    w <- log(rgamma(n, shape=k))
    y <- w / shape  + log(scale)
    exp(y) 
}

check.gengamma.orig <- function(shape, scale, k){
    ret <- TRUE
    if (missing(shape)) stop("shape parameter \"shape\" not given")
    if (missing(k)) stop("shape parameter \"k\" not given")
    if (any(shape <= 0)) {warning("Non-positive shape parameter \"shape\""); ret <- FALSE}
    if (any(scale <= 0)) {warning("Non-positive scale parameter"); ret <- FALSE}
    if (any(k <= 0)) {warning("Non-positive shape parameter \"k\""); ret <- FALSE}
    ret
}
