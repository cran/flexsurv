library(flexsurv)

## for local use
if (0) {
    library(survival)
    library(mvtnorm)
    library(muhaz)
    for (i in list.files("~/work/flexsurv/flexsurv/R", "*.R$"))
        source(paste("~/work/flexsurv/flexsurv/R/",i,sep=""))
}

test <- function(x, y, tol=1e-06) {
    stopifnot(isTRUE(all.equal(x, y, tol=tol)))
}

### TESTS WITH OVARIAN CANCER DATA FROM survival PACKAGE

## Basic GF fit -- doesn't converge -- "p" par not identifiable
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf",
                  control=list(trace=1,REPORT=1,maxit=10000,ndeps=rep(1e-06,4))))
## Basic GG fit
fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma")
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gengamma")
## GF with "p" fixed at 0
fitffix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf",
                       fixedpars=4, inits=c(NA,NA,NA,1e-05))
test(fitffix$res[1:3,"est"], fitg$res[1:3,"est"], tol=1e-03)
test(fitffix$res[1:3,2:3], fitg$res[1:3,2:3], tol=1e-02)

wt <- rep(1, nrow(ovarian)) # ; wt[c(1,3,5,7,9)] <- 10

## Weibull
fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull", weights=wt)
## Weibull with library(survival)
fitws <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull", weights=wt)
test(fitw$loglik, fitws$loglik[1], tol=1e-04)
test(fitws$scale, 1 / fitw$res["shape","est"], tol=1e-03)
test(as.numeric(coef(fitws)[1]), log(fitw$res["scale","est"]), tol=1e-03)
## Log-normal
fitln <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1,
                     data = ovarian, dist="lnorm")
## Gompertz
fitgo <- flexsurvreg(formula = Surv(futime, fustat) ~ 1,
                     data = ovarian, dist="gompertz",
                     control=list(trace=1,REPORT=1,reltol=1e-16))
fitgo
if (interactive()) {
    plot.flexsurvreg(fitf)
    lines.flexsurvreg(fitg, col="blue", lty.fit=2)
    lines.flexsurvreg(fitw, col="green", lty.fit=2)
    lines.flexsurvreg(fitln, col="purple", lty.fit=2)
    lines.flexsurvreg(fitgo, col="brown", lty.fit=2)
}

## Test distributions reducing to others with fixed pars
fitffix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf",
                       fixedpars=TRUE, inits=c(0,1,0,1))
## GG = GF with p -> 0
fitffix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf",
                       fixedpars=TRUE, inits=c(0,1,0,1e-08))
fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                       fixedpars=TRUE, inits=c(0,1,0))
test(fitgfix$loglik, fitffix$loglik, tol=1e-02)
## Weib = GG with q=1
fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                       fixedpars=TRUE, inits=c(6,0.8,1))
fitwfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull",
                       fixedpars=TRUE, inits=c(1/0.8,exp(6)))
test(fitwfix$loglik, fitgfix$loglik)
## Gamma = GG with sig=q
fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                       fixedpars=TRUE, inits=c(6,0.5,0.5))
fitgafix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gamma",
                       fixedpars=TRUE, inits=c(1/0.5^2,exp(-6)/0.5^2))
test(fitgafix$loglik,fitgfix$loglik)
## Log-normal = GG with q=0
fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                       fixedpars=TRUE, inits=c(6,0.8,0))
fitlfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",
                       fixedpars=TRUE, inits=c(6,0.8))
test(fitlfix$loglik,fitgfix$loglik)
## Compare with weib/lnorm fit from survreg
fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull",
                    inits=c(1/0.8,exp(6)))
fitw2 <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull")
test(1 / fitw2$scale, fitw$res["shape","est"], tol=1e-03)
test(as.numeric(coef(fitw2)[1]), log(fitw$res["scale","est"]), tol=1e-03)

## try a few optimisations with fixed pars and distributions reducing to others
fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma.orig",
                       fixedpars=3, inits=c(NA,NA,1))
fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull")
fitgfix
fitw
fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma.orig",
                       fixedpars=1, inits=c(1,NA,NA))
fitga <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gamma")
fitgfix
fitga
1/fitga$res["rate",]
fite <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="exp")
fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull", fixedpars=1)
fitw
1 / fite$res["rate",]

## Custom distribution
## Log-logistic
if (is.element("eha", installed.packages()[,1])) {
    library(eha)
    custom.llogis <- list(name="llogis",
                          pars=c("shape","scale"),
                          location="scale",
                          transforms=c(log, log),
                          inv.transforms=c(exp, exp),
                          inits=function(t){ c(1, median(t)) })
    fitll <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.llogis)
    fitll
    if (interactive()) {
        lines.flexsurvreg(fitll, col="pink", lty.fit=2)
    }
}

## SIMULATION TESTS. Simulate data from r and refit model to test
set.seed(12082012)
if (0) {
    sim <- rgenf(3000, 1.5, 1, -0.4, 0.6)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="genf", control=list(trace=1,REPORT=1))
    fit$res # OK

    sim <- rgengamma(3000, 1.5, 1, -0.4)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="gengamma", control=list(trace=1,REPORT=1))
    fit$res # OK

    sim <- rgenf.orig(3000, 1.5, 1, 0.4, 0.6)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="genf.orig", control=list(trace=1,REPORT=1,maxit=10000))
    fit$res # OK

    sim <- rgengamma.orig(3000, 1.5, 1, 0.4)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="gengamma.orig", control=list(trace=1,REPORT=1))
    fit$res # OK

    xg <- rgompertz(1000, 0.12, 4); hist(xg)
    flexsurvreg(Surv(xg, rep(1,1000)) ~ 1, dist="gompertz") ## OK. robust to starting values

    if (is.element("eha", installed.packages()[,1])) {
        library(eha)
        foo <- phreg(Surv(xg, rep(1,1000)) ~ 1, dist="gompertz") ## OK - names of parameters other way round, see dgompertz help
        xl <- rllogis(1000, 5.4, 0.1); hist(xl)
        flexsurvreg(Surv(xl, rep(1,1000)) ~ 1, dist=custom.llogis) ## OK. robust to starting values
    }
}

### TESTS WITH COVARIATES

fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(rx), data = ovarian, dist="weibull",
                  method="BFGS", control=list(trace=1,REPORT=1,maxit=10000,ndeps=rep(1e-06,3)))
fitg
if (interactive()) {
    plot.flexsurvreg(fitg, ci=TRUE)
    plot.flexsurvreg(fitg, X=rbind(c(0), c(1)), ci=TRUE, col="red")
    lines.flexsurvreg(fitg, X=rbind(c(1.1), c(1.2)), ci=TRUE, col="blue")
    plot.flexsurvreg(fitg, type="hazard")
    plot.flexsurvreg(fitg, type="cumhaz")
}

x <- rnorm(500,0,1)
sim <- rgenf(500, 1.5 - 0.2*x, 1, -0.4, 0.6)
dead <- as.numeric(sim<=30)
simt <- ifelse(sim<=30, sim, 30)
fit <- flexsurvreg(Surv(simt, dead) ~ x, dist="genf", control=list(trace=1,REPORT=1,maxit=10000))
fit
if (interactive()) {
    plot.flexsurvreg(fit)
    lines.flexsurvreg(fit, X=matrix(c(1,2),nrow=2))
    plot(fit)
    fit
    plot(fit, type="hazard", min.time=0, max.time=25)
    lines(fit, type="hazard", X=matrix(c(1,2),nrow=2))
    x2 <- factor(rbinom(500, 1, 0.5))
    fit <- flexsurvreg(Surv(simt, dead) ~ x + x2, dist="genf", control=list(trace=1,REPORT=1,maxit=10000))
    plot(fit)
    plot(fit, type="cumhaz")
    plot(fit, type="hazard", min.time=0, max.time=25)
    x3 <- factor(rbinom(500, 1, 0.5))
    fit <- flexsurvreg(Surv(simt, dead) ~ x2 + x3, dist="genf", control=list(trace=1,REPORT=1,maxit=10000))
    fit <- flexsurvreg(Surv(simt, dead) ~ x2, dist="genf", control=list(trace=1,REPORT=1,maxit=10000))
    plot(fit)
}

x2 <- factor(rbinom(500, 1, 0.5))
x3 <- rnorm(500,0,1)
sim <- rgengamma(500, 1.5 + 2*x3, 1, -0.4)
dead <- as.numeric(sim<=30)
simt <- ifelse(sim<=30, sim, 30)
fit <- flexsurvreg(Surv(simt, dead) ~ x3, dist="gengamma", control=list(trace=1,REPORT=1,maxit=10000))
fit <- flexsurvreg(Surv(simt, dead) ~ x + x2 + x3, dist="gengamma", control=list(trace=1,REPORT=1,maxit=10000))
fit <- flexsurvreg(Surv(simt, dead)[1:100] ~ x[1:100] + x2[1:100], dist="gengamma", control=list(trace=1,REPORT=1,maxit=10000), method="BFGS")
fit <- flexsurvreg(Surv(simt, dead)[1:100] ~ x[1:100], dist="gengamma", control=list(trace=1,REPORT=1,maxit=10000))
fit

## Covariates on auxiliary parameters
set.seed(11082012)
sim <- rgengamma(500, 1, exp(0.5 + 0.1*x3), -0.4)
dead <- as.numeric(sim<=30)
simt <- ifelse(sim<=30, sim, 30)
fit <- flexsurvreg(Surv(simt, dead) ~ sigma(x3), dist="gengamma", control=list(trace=1,REPORT=1,maxit=10000))
cl <- confint(fit)
stopifnot(cl[,1] < c(1, 0.5, -0.4, 0.1)  &  cl[,2] > c(1, 0.5, -0.4, 0.1) )


## Errors
if (0) {
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = c(1,2,3)))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = "foo"))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = c(1,2,3,-1)))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = c(1,-2,3,-1)))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", fixedpars = c(3,4,5,6,7)))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", fixedpars = "foo"))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl=-1))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl=1.1))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl=c(1,2)))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl="foo"))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="foo"))
try(fitf <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian))
}




### OTHER DATASETS (local developer use only)

if (0) {
    load("~/work/oral/tcrb.rda")
    ## Without covs
    fit <- flexsurvreg(Surv(survtime, dead2) ~ 1, data=tcrb, dist="genf")
    fit.gg <- flexsurvreg(Surv(survtime, dead2) ~ 1, data=tcrb, dist="gengamma")
    fit.sp2 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=2, data=tcrb, control=list(maxit=10000))
    fit.sp3 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=3, data=tcrb, control=list(maxit=10000))
    fit.sp4 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=4, data=tcrb, control=list(maxit=10000))
    fit.sp5 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=5, data=tcrb, method="BFGS", control=list(trace=1,REPORT=1,maxit=10000))
    ## min AIC with 4 knots, better than GF
    save(fit, fit.gg, fit.sp2, fit.sp3, fit.sp4, fit.sp5, file="~/work/flexsurv/tests/tcr.rda")
    plot(fit, ci=FALSE)
    lines(fit.gg, col="blue", ci=FALSE)
    lines(fit.sp4, col="green", ci=FALSE)

    ## With covs -- can't get arbitrarily good fit unless model 3 way interaction
    fitc.f <- flexsurvreg(Surv(survtime, dead2) ~ age_10 + sex + stage, data=tcrb, dist="genf")
    fitc.g <- flexsurvreg(Surv(survtime, dead2) ~ age_10 + sex + stage, data=tcrb, dist="gengamma")
    fitc.sp2 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=2, data=tcrb) # doesnt fit as well as gf/gg.  PH assumption vs AFT?
    fitc.sp3 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=3, data=tcrb)
    fitc.sp4 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=4, data=tcrb)
    fitc.spo2 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=2, scale="odds", data=tcrb) # PO doesn't fit better
    fitc.spn2 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=2, scale="normal", data=tcrb) # normal model is worse
    save(fitc.f, fitc.g, fitc.sp2, fitc.sp3, fitc.sp4, fitc.spo2, fitc.spn2, file="../../tests/tcrcov.rda")
    ## GF fits best, AFT assumption must be better than PH.

    plot(survfit(Surv(survtime, dead2) ~ stage, data=tcrb, subset=(tcrb$age_10=="50-59" & tcrb$sex=="female")))
    lines(fitc.f, X=rbind(c(1,0,0,0,1,0,0,0),
                  c(1,0,0,0,1,1,0,0),
                  c(1,0,0,0,1,0,1,0),
                  c(1,0,0,0,1,0,0,1)))

    fitc.f <- flexsurvreg(Surv(survtime, dead2) ~ stage, data=tcrb, dist="genf")
    plot(fitc.f)
    plot(fitc.f, type="hazard", min.time=0, max.time=10)
    plot(survfit(Surv(survtime, dead2) ~ stage, data=tcrb))
    lines(fitc.f, X=rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1)))
}

### calling flexsurvreg from within a function
### environment bug in early versions

f <- function(){
  ovarian2 <- ovarian
  fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian2, dist="gengamma")
  print(fitg)
  if(interactive()) plot(fitg)
  fitg <- flexsurvreg(formula = Surv(ovarian2$futime, ovarian2$fustat) ~ factor(ovarian2$rx), dist="gengamma",method="Nelder-Mead")
  if(interactive()) print(fitg)
  plot(fitg, ci=TRUE)
}
f()

if (interactive()) {
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma")
    plot(fitg)
    plot(fitg, type="cumhaz")
                                        # plot(fitg, type="hazard", min.time=0, max.time=1000)
    ## note can't change the ylim - need to use muhaz manually.
    fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data=ovarian, dist="gengamma")
    plot(fitg)

    ## does it work with lexical scoping
    f <- function(){
        ovarian2 <- ovarian
        g <- function(){
            fitw <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian2, dist="weibull")
            print(fitw)
            if(interactive()) plot(fitw)
            fitw <- flexsurvreg(formula = Surv(ovarian2$futime, ovarian2$fustat) ~ factor(ovarian2$rx), dist="weibull")
            if(interactive()) print(fitw)
            plot(fitw, ci=TRUE)
            plot(fitw, type="hazard", ci=TRUE)
        }
        g()
    }
    f()
}
## Left-truncation.
## Time passed as arg to initial values is stop - start,
## since, e.g. mean of trunc exponential dist is 1/lam + b, mean par plus trunc point
## time at risk in returned object is currently sum of (stop - start)
## default knot choice for spline - start + quantiles of log dt

if(0){
set.seed(12082012)
sim <- rgenf(3000, 1.5, 1, -0.4, 0.6)
dead <- as.numeric(sim<=30)
simt <- ifelse(sim<=30, sim, 30)
obs <- simt>3; simt <- simt[obs]; dead <- dead[obs]
fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="gengamma")
plot(fit, ci=FALSE, xlim=c(0,10))
fit <- flexsurvreg(Surv(rep(3, length(simt)), simt, dead) ~ 1, dist="gengamma")
lines(fit, ci=FALSE, col="blue") # truncated model fits truncated data better.
}
