### R code from vignette source 'flexsurv-examples.Rnw'

###################################################
### code chunk number 1: flexsurv-examples.Rnw:39-59
###################################################
library(flexsurv)

hgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
    dummy * hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}

HgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
    dummy * Hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}

custom.gengammaPH <- list(name="gengammaPH", 
                          pars=c("dummy","mu","sigma","Q"), location="dummy",
                          transforms=c(log, identity, log, identity),
                          inv.transforms=c(exp, identity, exp, identity),
                          inits=function(t){
                              lt <- log(t[t>0])
                              c(1, mean(lt), sd(lt), 0)
                          })
fs7 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, 
                   dist=custom.gengammaPH, fixedpars=1)


###################################################
### code chunk number 2: flexsurv-examples.Rnw:84-127
###################################################

fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ group + sigma(group), 
                   data=bc, dist="gengamma")
B <- 5000
t <- seq(0.1, 8, by=0.1)

hrAFT.est <-
    summary.flexsurvreg(fs2, t=t, type="hazard",
                        newdata=data.frame(group="Medium"),ci=FALSE)[[1]][,"est"] /
    summary.flexsurvreg(fs2, t=t, type="hazard",
                        newdata=data.frame(group="Good"),ci=FALSE)[[1]][,"est"]
pars <- normboot.flexsurvreg(fs2, B=B, newdata=data.frame(group=c("Good","Medium")))
hrAFT <- matrix(nrow=B, ncol=length(t))
for (i in seq_along(t)){ 
    haz.medium.rep <- do.call(hgengamma, c(list(t[i]), as.data.frame(pars[[2]])))
    haz.good.rep <- do.call(hgengamma, c(list(t[i]), as.data.frame(pars[[1]])))
    hrAFT[,i] <- haz.medium.rep / haz.good.rep
}
hrAFT <- apply(hrAFT, 2, quantile, c(0.025, 0.975))

hrPH.est <-
    summary.flexsurvreg(fs7, t=t, type="hazard",
                        newdata=data.frame(group="Medium"),ci=FALSE)[[1]][,"est"] /
    summary.flexsurvreg(fs7, t=t, type="hazard",
                        newdata=data.frame(group="Good"),ci=FALSE)[[1]][,"est"]
pars <- normboot.flexsurvreg(fs7, B=B, newdata=data.frame(group=c("Good","Medium")))
hrPH <- matrix(nrow=B, ncol=length(t))
for (i in seq_along(t)){ 
    haz.medium.rep <- do.call(hgengammaPH, c(list(t[i]), as.data.frame(pars[[2]])))
    haz.good.rep <- do.call(hgengammaPH, c(list(t[i]), as.data.frame(pars[[1]])))
    hrPH[,i] <- haz.medium.rep / haz.good.rep
}
hrPH <- apply(hrPH, 2, quantile, c(0.025, 0.975))

plot(t, hrAFT[1,], type="l", ylim=c(0, 10), col="red", xlab="Years",
     ylab="Hazard ratio (Medium / Good)", lwd=1, lty=2)
lines(t, hrAFT[2,], col="red", lwd=1, lty=2)
lines(t, hrPH[1,], col="darkgray", lwd=1, lty=2)
lines(t, hrPH[2,], col="darkgray", lwd=1, lty=2)
lines(t, hrAFT.est, col="red", lwd=2)
lines(t, hrPH.est, col="darkgray", lwd=2)
legend("topright", lwd=c(2,2), col=c("red","darkgray"), bty="n",
       c("Generalized gamma: standard AFT", "Generalized gamma: proportional hazards"))


