### R code from vignette source 'flexsurv.Rnw'

###################################################
### code chunk number 1: flexsurv.Rnw:186-188
###################################################
library(flexsurv)
bc[1:2,]


###################################################
### code chunk number 2: flexsurv.Rnw:194-195
###################################################
fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")


###################################################
### code chunk number 3: flexsurv.Rnw:213-214
###################################################
fs1


###################################################
### code chunk number 4: flexsurv.Rnw:220-221
###################################################
survreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")


###################################################
### code chunk number 5: flexsurv.Rnw:319-322
###################################################
fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="gengamma")
fs3 <- flexsurvreg(Surv(recyrs, censrec) ~ group + sigma(group), 
                   data=bc, dist="gengamma")


###################################################
### code chunk number 6: surv
###################################################
plot(fs1, col="gray", lwd.obs=2, xlab="Years", ylab="Recurrence-free survival")
lines(fs2, col="red", lty=2)
lines(fs3, col="red")
legend("bottomleft", col=c("black","gray","red","red"), 
       lty=c(1,1,2,1), bty="n", lwd=rep(2,4),
       c("Kaplan-Meier","Weibull","Generalized gamma (AFT)",
         "Generalized gamma (time-varying)"))


###################################################
### code chunk number 7: haz
###################################################
plot(fs1, type="hazard", col="gray", lwd.obs=2, xlab="Years", ylab="Hazard")
lines(fs2, type="hazard", col="red", lty=2)
lines(fs3, type="hazard", col="red")
legend("topright", col=c("black","gray","red","red"),
       lty=c(1,1,2,1),  bty="n", lwd=rep(2,4),
       c("Kernel density estimate","Weibull","Gen. gamma (AFT)",
         "Gen. gamma (time-varying)"))


###################################################
### code chunk number 8: flexsurv.Rnw:450-454
###################################################
median.weibull <- function(shape, scale) { 
    qweibull(0.5, shape=shape, scale=scale) 
}
summary(fs1, fn=median.weibull, t=1, B=10000)


###################################################
### code chunk number 9: flexsurv.Rnw:588-593
###################################################
library(eha)
custom.llogis <- list(name="llogis",  pars=c("shape","scale"), location="scale",
                      transforms=c(log, log), inv.transforms=c(exp, exp),
                      inits=function(t){ c(1, median(t)) })
fs4 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.llogis)


###################################################
### code chunk number 10: flexsurv.Rnw:611-617
###################################################
dmakeham3 <- function(x, shape1, shape2, scale, ...)  {
    dmakeham(x, shape=c(shape1, shape2), scale=scale, ...)
}
pmakeham3 <- function(q, shape1, shape2, scale, ...)  {
    pmakeham(q, shape=c(shape1, shape2), scale=scale, ...)
}


###################################################
### code chunk number 11: flexsurv.Rnw:624-626
###################################################
dmakeham3 <- Vectorize(dmakeham3) 
pmakeham3 <- Vectorize(pmakeham3)


###################################################
### code chunk number 12: flexsurv.Rnw:629-630
###################################################
pmakeham3(c(0, 1, 1, Inf), 1, c(1, 1, 2, 1), 1)


###################################################
### code chunk number 13: flexsurv.Rnw:646-662
###################################################
detach("package:eha")
hweibullPH <- function(x, shape, scale = 1, log=FALSE){
    hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}
HweibullPH <- function(x, shape, scale=1, log=FALSE){
    Hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}
custom.weibullPH <- list(name="weibullPH", 
                         pars=c("shape","scale"), location="scale",
                         transforms=c(log, log), inv.transforms=c(exp, exp),
                         inits = function(t){
                             c(1, median(t[t>0]) / log(2))
                         })
fs6 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.weibullPH)
fs6$res["scale","est"] ^ {-1/fs6$res["shape","est"]}
- fs6$res["groupMedium","est"] / fs6$res["shape","est"]


###################################################
### code chunk number 14: flexsurv.Rnw:781-783
###################################################
sp1 <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, 
                      scale="odds")


###################################################
### code chunk number 15: flexsurv.Rnw:787-789
###################################################
sp2 <- flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(group),
                      data=bc, k=1, scale="odds")


###################################################
### code chunk number 16: splinehaz
###################################################
plot(sp1, type="hazard", ylim=c(0, 0.5), xlab="Years", ylab="Hazard")
lines(sp2, type="hazard", col="red", lty=2)
lines(fs2, type="hazard", col="blue")
legend("topright", col=c("black","red","red","blue"), lty=c(1,1,2,1), lwd=rep(2,4),
       c("Kernel density estimate","Spline (proportional odds)",
         "Spline (time-varying)","Generalized gamma (time-varying)"))


###################################################
### code chunk number 17: flexsurv.Rnw:812-816
###################################################
sp3 <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="hazard")
sp3$res[c("groupMedium","groupPoor"),c("est","se")]
cox3 <- coxph(Surv(recyrs, censrec) ~ group, data=bc)
coef(summary(cox3))[,c("coef","se(coef)")]


###################################################
### code chunk number 18: flexsurv.Rnw:820-826
###################################################
res <- t(sapply(list(fs1, fs2, fs3, fs4, sp1, sp2), 
                function(x)rbind(-2*x$loglik, x$npars, x$AIC)))
rownames(res) <- c("Weibull (fs1)","Generalized gamma (fs2)",
                   "Generalized gamma (fs3)","Log-logistic (fs4)",
                   "Spline (sp1)", "Spline (sp2)")
colnames(res) <- c("-2 log likelihood","Parameters","AIC")


###################################################
### code chunk number 19: flexsurv.Rnw:829-830
###################################################
res


###################################################
### code chunk number 20: flexsurv.Rnw:869-871
###################################################
gamma <- sp1$res[c("gamma0","gamma1","gamma2"),"est"]
1 - psurvspline(5, gamma=gamma, knots=sp1$knots)


###################################################
### code chunk number 21: flexsurv.Rnw:877-879
###################################################
pfn <- unroll.function(psurvspline, gamma=0:2)
1 - pfn(5, gamma0=gamma[1], gamma1=gamma[2], gamma2=gamma[3], knots=sp1$knots)


###################################################
### code chunk number 22: flexsurv.Rnw:913-922
###################################################
hsurvspline.lh <- function(x, gamma, knots){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
    lg <- nrow(gamma) # return vector of length of longest argument
    nret <- max(length(x), lg)
    gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
    x <- rep(x, length=nret)
    loghaz <- rowSums(basis(knots, log(x)) * gamma)
    exp(loghaz)
}


###################################################
### code chunk number 23: flexsurv.Rnw:928-929
###################################################
hsurvspline.lh3 <- unroll.function(hsurvspline.lh, gamma=0:2)


###################################################
### code chunk number 24: flexsurv.Rnw:931-937
###################################################
custom.hsurvspline.lh3 <- list(
    name = "survspline.lh3",
    pars = c("gamma0","gamma1","gamma2"),
    location = c("gamma0"),
    transforms = rep(c(identity), 3), inv.transforms=rep(c(identity), 3)
    )


###################################################
### code chunk number 25: flexsurv.Rnw:942-944
###################################################
dtime <- log(bc$recyrs)[bc$censrec==1]
ak <- list(knots=quantile(dtime, c(0, 0.5, 1)))


###################################################
### code chunk number 26: flexsurv.Rnw:954-958 (eval = FALSE)
###################################################
## sp4 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, aux=ak,
##                    inits=c(0, 0, 0, 0, 0), dist=custom.hsurvspline.lh3, 
##                    method="L-BFGS-B", lower=c(-Inf,-Inf,-0.5), upper=c(Inf,Inf,0.5),
##                    control=list(trace=1,REPORT=1))


###################################################
### code chunk number 27: flexsurv.Rnw:1103-1104
###################################################
bosms3[18:22,]


###################################################
### code chunk number 28: flexsurv.Rnw:1146-1152
###################################################
crexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, 
                     dist="exp")
crwei <- flexsurvreg(Surv(years, status) ~ trans + shape(trans), 
                     data=bosms3, dist="weibull")
cfwei <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + shape(trans), 
                     data=bosms3, dist="weibull")


###################################################
### code chunk number 29: flexsurv.Rnw:1162-1164
###################################################
crcox <- coxph(Surv(years, status) ~ strata(trans), data=bosms3)
cfcox <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=bosms3)


###################################################
### code chunk number 30: flexsurv.Rnw:1199-1203
###################################################
library(mstate)
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
mrcox <- msfit(crcox, trans=tmat)
mfcox <- msfit(cfcox, trans=tmat)


###################################################
### code chunk number 31: flexsurv.Rnw:1223-1227
###################################################
tgrid <- seq(0,14,by=0.1)
mrwei <- msfit.flexsurvreg(crwei, t=tgrid, trans=tmat)
mrexp <- msfit.flexsurvreg(crexp, t=tgrid, trans=tmat)
mfwei <- msfit.flexsurvreg(cfwei, t=tgrid, trans=tmat)


###################################################
### code chunk number 32: cumhaz
###################################################
cols <- c("black","red","blue")
plot(mrcox, xlab="Years after baseline", lwd=3, xlim=c(0,14), cols=cols)
for (i in 1:3){
    lines(tgrid, mrexp$Haz$Haz[mrexp$Haz$trans==i], col=cols[i], lty=2, lwd=2)
    lines(tgrid, mrwei$Haz$Haz[mrwei$Haz$trans==i], col=cols[i], lty=3, lwd=2)
}
lines(mfcox$Haz$time[mfcox$Haz$trans==3], mfcox$Haz$Haz[mfcox$Haz$trans==3],
      type="s", col="darkgreen", lty=1, lwd=2)
lines(tgrid, mfwei$Haz$Haz[mfwei$Haz$trans==3], col="darkgreen", lty=3, lwd=2)
legend("topleft", inset=c(0,0.2), lwd=2, col=c("darkgreen"), 
       c("2 -> 3 (clock-forward)"), bty="n")
legend("topleft", inset=c(0,0.3), c("Non-parametric","Exponential","Weibull"),
       lty=c(1,2,3), lwd=c(3,2,2), bty="n")


###################################################
### code chunk number 33: flexsurv.Rnw:1277-1278
###################################################
pmatrix.fs(cfwei, t=c(5,10), trans=tmat)


###################################################
### code chunk number 34: flexsurv.Rnw:1305-1307 (eval = FALSE)
###################################################
## pmatrix.simfs(crwei, trans=tmat, t=5)
## pmatrix.simfs(crwei, trans=tmat, t=10)


###################################################
### code chunk number 35: flexsurv.Rnw:1331-1333
###################################################
ptc <- probtrans(mfcox, predt=0, direction="forward")[[1]]
ptc[c(165, 193),]


###################################################
### code chunk number 36: flexsurv.Rnw:1341-1343
###################################################
ptw <- probtrans(mfwei, predt=0, direction="forward")[[1]]
ptw[ptw$time %in% c(5,10),]


###################################################
### code chunk number 37: flexsurv.Rnw:1365-1367 (eval = FALSE)
###################################################
## mssample(mrcox$Haz, trans=tmat, clock="reset", M=1000, tvec=c(5, 10))
## mssample(mrwei$Haz, trans=tmat, clock="reset", M=1000, tvec=c(5, 10))


