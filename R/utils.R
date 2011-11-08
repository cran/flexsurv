### Quantile function "qdist" for a generic distribution "dist".
### Works using an interval search for the solution of pdist(q) = p.
### Requires a probability function "pdist" for the same distribution
### in the working environment.

qgeneric <- function(pdist, p, lower.tail = TRUE, log.p = FALSE, ...)
{
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    ret <- numeric(length(p))
    ret[p == 1] <- Inf
    ret[p == 0] <- -Inf
    ret[p < 0 | p > 1] <- NaN
    ind <- (p > 0 & p < 1)
    if (any(ind)) {
        hind <- seq(along=p)[ind]
        h <- function(y) {
            args <- list(...)
            args$q <- y
            (do.call(pdist, args) - p)[hind[i]]
        }
        ptmp <- numeric(length(p[ind]))
        for (i in 1:length(p[ind])) {
            interval <- c(-1, 1)
            while (h(interval[1])*h(interval[2]) >= 0) { 
              interval <- interval + c(-1,1)*0.5*(interval[2]-interval[1])
          }
            ptmp[i] <- uniroot(h, interval, tol=.Machine$double.eps)$root
        }
        ret[ind] <- ptmp
    }
    if (any(is.nan(ret))) warning("NaNs produced")
    ret
}
