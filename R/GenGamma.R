## Log-gamma or generalized gamma distribution (parameterisation as in Farewell and Prentice, Technometrics 1977)

dgengamma <- function(x, mu=0, sigma=1, Q, log=FALSE) {
    n <- max(length(x),length(mu),length(sigma),length(Q))
    x <- rep(x, length=n)
    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    Q <- rep(Q, length=n)
    ret <- numeric(n)
    ret[!check.gengamma(mu=mu, sigma=sigma, Q=Q)] <- NaN
    if (all(is.nan(ret))) return(ret);
    ret[!is.nan(ret) & (x<=0)] <- if (log) -Inf else 0
    ind <- !is.nan(ret) & (x>0)
    mu <- mu[ind]; sigma <- sigma[ind]; Q <- Q[ind]; xx <- x[ind]
    logdens <- numeric(length(xx))
    logdens[Q==0] <- dlnorm(xx[Q==0], mu[Q==0], sigma[Q==0], log=TRUE)
    qn0 <- Q!=0
    if (any(qn0)) {
        xx <- xx[qn0]; mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
        y <- log(xx)
        w <- ((y - mu)/sigma)
        logdens[qn0] <- -log(sigma*xx) + log(abs(Q)) + (Q^-2)*log(Q^-2) + Q^-2*(Q*w - exp(Q*w)) - lgamma(Q^-2)
    }
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}

pgengamma <- function(q, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE) {
    n <- max(length(q),length(mu),length(sigma),length(Q))
    q <- rep(q, length=n)
    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    Q <- rep(Q, length=n)
    ret <- numeric(n)
    ret[!check.gengamma(mu=mu, sigma=sigma, Q=Q)] <- NaN
    if (all(is.nan(ret))) return(ret);
    q[q<0] <- 0
    ind <- !is.nan(ret)
    q <- q[ind]; mu <- mu[ind]; sigma <- sigma[ind]; Q <- Q[ind]
    prob <- numeric(length(q))
    prob[Q==0] <- plnorm(q[Q==0], mu[Q==0], sigma[Q==0])
    qn0 <- Q!=0
    if (any(qn0)) {
        q <- q[qn0]; mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
        y <- log(q)
        w <- ((y - mu)/sigma)
        expnu <- exp(Q*w)*Q^-2
        prob[qn0] <- ifelse(Q > 0, pgamma(expnu, Q^-2), 1 - pgamma(expnu, Q^-2))
    }
    if (!lower.tail) prob <- 1 - prob
    if (log.p) prob <- log(prob)
    ret[ind] <- prob
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
    n <- max(length(p),length(mu),length(sigma),length(Q))
    p <- rep(p, length=n)
    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    Q <- rep(Q, length=n)
    ret <- numeric(n)
    ret[!check.gengamma(mu=mu, sigma=sigma, Q=Q)] <- NaN
    if (all(is.nan(ret))) return(ret);
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    ind <- !is.nan(ret)
    p <- p[ind]; mu <- mu[ind]; sigma <- sigma[ind]; Q <- Q[ind]
    p[Q<0] <- 1 - p[Q<0]
    ret[ind] <- numeric(sum(ind))
    ret[ind][Q==0] <- qlnorm(p[Q==0], mu[Q==0], 1/sigma[Q==0]^2)
    qn0 <- Q!=0
    p <- p[qn0]; mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
    ret[ind][qn0] <- exp(mu + sigma*(log(Q^2*qgamma(p, 1/Q^2, 1)) / Q))
    ret
}

rgengamma <- function(n, mu=0, sigma=1, Q) {
    if (length(n) > 1) n <- length(n)
    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    Q <- rep(Q, length=n)
    ret <- numeric(n)
    ret[!check.gengamma(mu=mu, sigma=sigma, Q=Q)] <- NaN
    if (all(is.nan(ret))) return(ret);
    ind <- !is.nan(ret)
    mu <- mu[ind]; sigma <- sigma[ind]; Q <- Q[ind]
    ret[ind][Q==0] <- rlnorm(n, mu, 1/sigma^2)
    qn0 <- Q!=0
    if (any(qn0)) {
        mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
        w <- log(Q^2*rgamma(n, 1/Q^2, 1)) / Q
        ret[ind][qn0] <- exp(mu + sigma*w)
    }
    ret
}

check.gengamma <- function(mu, sigma, Q){
    ret <- rep(TRUE, length(mu))
    if (missing(Q)) stop("shape parameter \"Q\" not given")
    if (any(sigma <= 0)) {warning("Non-positive scale parameter \"sigma\""); ret[sigma<=0] <- FALSE}
    ret
}

dgengamma.orig <- function(x, shape, scale=1, k, log=FALSE){
    n <- max(length(x),length(shape),length(scale),length(k))
    x <- rep(x, length=n)
    shape <- rep(shape, length=n)
    scale <- rep(scale, length=n)
    k <- rep(k, length=n)
    ret <- numeric(n)
    ret[!check.gengamma.orig(shape=shape, scale=scale, k=k)] <- NaN
    if (all(is.nan(ret))) return(ret);
    ret[!is.nan(ret) & (x<=0)] <- if (log) -Inf else 0
    ind <- !is.nan(ret) & (x>0)
    x <- x[ind]; shape <- shape[ind]; scale <- scale[ind]; k <- k[ind]
    logdens <- log(shape) - lgamma(k) + (shape*k - 1)*log(x) - shape*k*log(scale) - (x/scale)^shape
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}

pgengamma.orig <- function(q, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE) {
    n <- max(length(q),length(shape),length(scale),length(k))
    q <- rep(q, length=n)
    shape <- rep(shape, length=n)
    scale <- rep(scale, length=n)
    k <- rep(k, length=n)
    ret <- numeric(n)
    ret[!check.gengamma.orig(shape=shape, scale=scale, k=k)] <- NaN
    if (all(is.nan(ret))) return(ret);
    q[q<0] <- 0
    ind <- !is.nan(ret)
    q <- q[ind]; shape <- shape[ind]; scale <- scale[ind]; k <- k[ind]
    y <- log(q)
    w <- (y - log(scale))*shape
    prob <- pgamma(exp(w), shape=k)
    if (!lower.tail) prob <- 1 - prob
    if (log.p) prob <- log(prob)
    ret[ind] <- prob
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
    n <- max(length(p),length(shape),length(scale),length(k))
    p <- rep(p, length=n)
    shape <- rep(shape, length=n)
    scale <- rep(scale, length=n)
    k <- rep(k, length=n)
    ret <- numeric(n)
    ret[!check.gengamma.orig(shape=shape, scale=scale, k=k)] <- NaN
    if (all(is.nan(ret))) return(ret);
    ind <- !is.nan(ret)
    p <- p[ind]; shape <- shape[ind]; scale <- scale[ind]; k <- k[ind]
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    w <- log(qgamma(p, shape=k))
    y <- w / shape  + log(scale)
    ret[ind] <- exp(y)
    ret
}

rgengamma.orig <- function(n, shape, scale=1, k) {
    if (length(n) > 1) n <- length(n)
    shape <- rep(shape, length=n)
    scale <- rep(scale, length=n)
    k <- rep(k, length=n)
    ret <- numeric(n)
    ret[!check.gengamma.orig(shape=shape, scale=scale, k=k)] <- NaN
    if (all(is.nan(ret))) return(ret);
    ind <- !is.nan(ret)
    shape <- shape[ind]; scale <- scale[ind]; k <- k[ind]
    w <- log(rgamma(n, shape=k))
    y <- w / shape  + log(scale)
    ret[ind] <- exp(y)
    ret
}

check.gengamma.orig <- function(shape, scale, k){
    ret <- rep(TRUE, length(shape))
    if (missing(shape)) stop("shape parameter \"shape\" not given")
    if (missing(k)) stop("shape parameter \"k\" not given")
    if (any(shape <= 0)) {warning("Non-positive shape parameter \"shape\""); ret[shape<=0] <- FALSE}
    if (any(scale <= 0)) {warning("Non-positive scale parameter"); ret[scale<=0] <- FALSE}
    if (any(k <= 0)) {warning("Non-positive shape parameter \"k\""); ret[k<=0] <- FALSE}
    ret
}
