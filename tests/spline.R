### SPLINES
library(flexsurv)
data(bc)

## for local use 
if (0) {
    library(survival)
    library(muhaz)
    for (i in list.files("../R", "*.R$")) 
        source(paste("../R/",i,sep=""))
    bc <- read.table("../data/bc.txt")
}

test <- function(x, y, tol=1e-06) {
    stopifnot(isTRUE(all.equal(x, y, tol=tol)))
}

bc$foo <- factor(sample(1:3, nrow(bc), replace=TRUE))
bc$cts <- rnorm(nrow(bc))
bc$recyrs <- bc$rectime/365

## Weibull
spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=0)
spl
spl$loglik  +  sum(log(bc$recyrs[bc$censrec==1])) # OK same as Stata streg , and stpm with hazard = TRUE 
-2*(spl$loglik  +  sum(log(bc$recyrs[bc$censrec==1])))
if (interactive()){ 
    plot(spl)
    plot(spl, type="cumhaz")
    plot(spl, type="haz")
}

## Best-fitting model for breast cancer example in paper 
## gamma -3.451(0.203), 2.915(0.298), 0.191(0.044) 
spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="odds")
spl
spl$loglik
spl$loglik  +   sum(log(bc$recyrs[bc$censrec==1]))
spl$AIC  +   sum(log(bc$recyrs[bc$censrec==1]))
if (interactive()){ 
    plot(spl)
    plot(spl, type="cumhaz")
    plot(spl, type="haz")
}


splh <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="hazard")
spln <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="normal")
if (interactive()){
    plot(spl)
    lines(splh, col="blue")
    lines(spln, col="green")
}


### Test reduction to weibull
### what are pars? log(H(t)) = g0 + g1 log(t) + bz 
### H(t) = exp(g0) t^g1 exp(bz) 
### S(t)  = exp(-H(t))   = exp( - exp(g0) t^g1 exp(bz) )
### pweibull has par exp( - (x/b) ^ a)
### a = g1, 1/b^a = exp(g0 + bz) = b ^ -a  = exp(-a log b)
### g0 + bz = - a (log b0 + bz)

wei <- survreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")
wei.base <- survreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull")
a <- 1/wei$scale
b1 <- exp(coef(wei)[1]); b2 <- exp(coef(wei)[1]+coef(wei)[2]); b3 <- exp(coef(wei)[1]+coef(wei)[3])
a.base <- 1/wei.base$scale
b.base <- exp(coef(wei.base[1]))

## Compare three implementations of the Weibull, with and without covs 
fit <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull", fixedpars=FALSE,
                   inits=c(a,b1,coef(wei)[2:3]))
fit$loglik
spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=0, 
                      inits=c(-a*log(b1), a, -a*coef(wei)[2:3]), fixedpars=FALSE)
spl$"loglik" 
test(fit$loglik, spl$loglik)
test(fit$loglik, wei$loglik[2])

fit <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull", fixedpars=FALSE,
                   inits=c(a,b1))
spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0, 
                      inits=c(log(1 / b.base^a.base), a.base), fixedpars=FALSE)
test(fit$loglik, spl$loglik)
test(fit$loglik, wei.base$loglik[1])


### Test log-logistic reduction 
if (is.element("eha", installed.packages()[,1])) { 
    library(eha)
    custom.llogis <- list(name="llogis",
                          pars=c("shape","scale"),
                          location="scale",
                          transforms=c(log, log),
                          inv.transforms=c(exp, exp),
                          inits=function(t){ c(1, median(t)) })
    fitll <- flexsurvreg(formula = Surv(recyrs, censrec) ~ 1, data = bc, dist=custom.llogis)
    fitll 
    fitsp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0, scale="odds")
    test(fitsp$loglik, fitll$loglik)
    test(1/fitll$res["scale",1]^fitll$res["shape",1], exp(fitsp$res["gamma0",1]), tol=1e-02)
    test(fitsp$res["gamma1",1], fitll$res["shape",1], tol=1e-02)
    if (interactive()) { 
        lines.flexsurvreg(fitll, col.fit="pink", lty.fit=2)
    }
}

### Test log-normal reduction 
fitln <- flexsurvreg(formula = Surv(recyrs, censrec) ~ 1, data = bc, dist="lnorm")
fitln
fitsp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0, scale="normal")
test(fitsp$res["gamma0",1], -fitln$res["meanlog",1]/fitln$res["sdlog",1], tol=1e-02)
test(fitsp$res["gamma1",1], 1 /fitln$res["sdlog",1], tol=1e-02)
