context("Multi-state modelling and prediction")

bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp") 
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
tgrid <- seq(0,14,by=0.1)
bwei <- flexsurvreg(Surv(years, status) ~ trans + shape(trans), data=bosms3, dist="weibull")
bspl <- flexsurvspline(Surv(years, status) ~ trans + gamma1(trans), data=bosms3, k=3)
bexp.markov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")
bln.markov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="lnorm")

test_that("msfit.flexsurvreg",{
    mexp <- msfit.flexsurvreg(bexp, t=0.01, trans=tmat, tvar="trans")
    summ <- summary(bexp, t=0.01, type="cumhaz", ci=FALSE, 
                    newdata=list(trans=factor(1:3, levels=1:3)))
    summ <- as.numeric(unlist(lapply(summ, function(x)x$est[x$time==0.01])))   
    expect_equal(mexp$Haz$Haz[mexp$Haz$time==0.01], summ)
    mwei <- msfit.flexsurvreg(bwei, t=c(0.01, 0.02), trans=tmat, tvar="trans", B=10)
    mspl <- msfit.flexsurvreg(bspl, t=c(0.01, 0.02), trans=tmat, tvar="trans", B=10)
})

## With covariates
suppressWarnings(RNGversion("3.5.0"))
set.seed(1)
bosms3$x <- rnorm(nrow(bosms3))
bexp.cov <- flexsurvreg(Surv(years, status) ~ trans + x, data=bosms3, dist="exp")
bexp.markov.cov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + x, data=bosms3, dist="exp")

test_that("newdata in msfit.flexsurvreg",{
    expect_error({
        msfit.flexsurvreg(bexp.cov, newdata=list(x=1), t=c(0,5,10), trans=tmat, variance=FALSE)
        msfit.flexsurvreg(bexp.cov, newdata=list(x=2), t=c(0,5,10), trans=tmat, variance=FALSE)
        msfit.flexsurvreg(bexp.cov, newdata=list(x=c(1,2,3)), t=c(0,5,10), trans=tmat, variance=FALSE)
    }, NA)
})

test_that("Errors in msfit.flexsurvreg",{
    expect_error(msfit.flexsurvreg(bexp.cov, t=seq(0,150,10), trans=tmat), "Value.* of covariate.* .+ not supplied")
    expect_error(msfit.flexsurvreg(bexp.cov, t=seq(0,150,10), trans=tmat, tvar="foo"), "variable .* not in model")
    expect_error(msfit.flexsurvreg(bexp.cov, newdata=list(x=c(1,2)), t=c(0,5,10), trans=tmat, variance=FALSE), "length of variables .+ must be")   
})

test_that("pmatrix.fs",{
    pmat <- pmatrix.fs(bexp.markov, t=c(5,10), trans=tmat)
    expect_equal(pmat$"5"[1,2], 0.267218506920585, tolerance=1e-04)
    pmat <- pmatrix.fs(bexp.markov.cov, t=c(5,10), trans=tmat, newdata=list(x=1))
    expect_equal(pmat$"5"[1,2], 0.259065437633427, tolerance=1e-04)
})

test_that("totlos.fs for flexsurvreg objects",{
    tl <- totlos.fs(bexp.markov, t=c(5), trans=tmat)
    expect_equal(as.numeric(tl), c(2.89231556324412, 0, 0, 1.06822543404334, 2.77639174263866, 0, 1.03945900271255, 2.22360825736133, 5), tolerance=1e-06)
    tl <- totlos.fs(bexp.markov.cov, t=c(5), trans=tmat, newdata=list(x=1))
    expect_equal(as.numeric(tl), c(2.76115934740386, 0, 0, 1.08844199049896, 2.64568022609759, 0, 1.15039866209718, 2.35431977390241, 5))
    tl <- totlos.fs(bexp.markov, t=c(5,10), trans=tmat)
    expect_equal(as.numeric(tl[[1]]), c(2.89231557124359, 0, 0, 1.06822540869797, 2.77639177178247, 0, 1.03945902005845, 2.22360822821753, 5))
    tl <- totlos.fs(bexp.markov.cov, t=c(5,10), trans=tmat, newdata=list(x=1))
    expect_equal(as.numeric(tl[[1]]),c(2.76115934740386, 0, 0, 1.08844199049896, 2.64568022609759, 0, 1.15039866209718, 2.35431977390241, 5))
    attr(tl, "P")
})

test_that("pmatrix.simfs for flexsurvreg objects",{
  expect_error({
    pmatrix.simfs(bwei, t=5, trans=tmat, M=100)
    pmatrix.simfs(bexp.cov, t=5, trans=tmat, newdata=list(x=1), M=100)
  }, NA)
  dimnames(tmat) <- rep(list(paste("State",1:3)),2)
  set.seed(1)
  p5 <- pmatrix.simfs(bexp, t=5, trans=tmat, M=1000)
  p10 <- pmatrix.simfs(bexp, t=10, trans=tmat, M=1000)
  pboth <- pmatrix.simfs(bexp, t=c(0, 5,10), trans=tmat)
  expect_equivalent(p5, pboth[,,"5"], tol=0.1)
  (p5 <- pmatrix.simfs(bexp, t=5, trans=tmat, ci=TRUE, M=100, B=3))
  expect_gt(p5[1,1], attr(p5,"lower")[1,1])
  set.seed(1)
  (pboth <- pmatrix.simfs(bexp, t=c(0, 5,10), trans=tmat, ci=TRUE, M=1000, B=3))
  expect_gt(pboth[1,1,"5"], attr(pboth,"lower")[1,1,"5"])
  set.seed(1)
  (pt <- pmatrix.simfs(bexp, t=c(0, 5,10), trans=tmat, ci=TRUE, M=1000, B=3, tidy=TRUE))
  expect_equal(pboth[1,2,"5"], pt[pt$from=="State 1" & pt$to=="State 2" & pt$t==5, "p"])
  expect_equal(attr(pboth,"lower")[1,2,"5"], 
               pt[pt$from=="State 1" & pt$to=="State 2" & pt$t==5, "lower"])
  expect_equal(attr(pboth,"upper")[1,2,"5"], 
               pt[pt$from=="State 1" & pt$to=="State 2" & pt$t==5, "upper"])
})

test_that("totlos.simfs for flexsurvreg objects",{
    expect_error({
        totlos.simfs(bexp, t=5, trans=tmat, M=100)
        totlos.simfs(bwei, t=5, trans=tmat, M=100)
        totlos.simfs(bexp.cov, t=5, trans=tmat, newdata=list(x=1), M=100)
        totlos.simfs(bexp, t=5, trans=tmat, M=100)
        totlos.simfs(bwei, t=5, trans=tmat, M=100)
        totlos.simfs(bexp.cov, t=5, trans=tmat, newdata=list(x=1), M=100)
    }, NA)
})

### List format for independent transition-specific models

bwei.list <- bweic.list <- bweim.list <- bexpc.list <- vector(3, mode="list")
for (i in 1:3) {
    bwei.list[[i]] <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==i),
                                   data = bosms3, dist = "weibull")
    bweic.list[[i]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==i),
                                   data = bosms3, dist = "weibull")
    bweim.list[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ 1, subset=(trans==i),
                               data=bosms3, dist="weibull")
    bexpc.list[[i]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==i), data=bosms3, dist="exp")
}

test_that("multistate output functions for list of models objects", {
    set.seed(1)
    totlos.simfs(bwei.list, t=5, trans=tmat, M=10)
    totlos.simfs(bweic.list, t=5, trans=tmat, M=100, newdata=list(x=0))
    totlos.simfs(bweic.list, t=5, trans=tmat, M=100, newdata=list(x=0), tcovs="x")
    
    pmatrix.simfs(bwei.list, t=5, trans=tmat, M=100)
    pmatrix.simfs(bweic.list, t=10, trans=tmat, M=10000, newdata=list(x=0), tcovs="x")
    
    set.seed(1)
    pnt <- pmatrix.simfs(bweic.list, t=5, trans=tmat, M=100, newdata=list(x=0))
    set.seed(1)
    pt <- pmatrix.simfs(bweic.list, t=c(5), trans=tmat, M=100, newdata=list(x=0),
                        tidy=TRUE)
    expect_equal(pnt[1,2], pt$p[pt$from=="1" & pt$to=="2" & pt$t==5])
    
    pmatrix.fs(bweim.list, t=5, trans=tmat)
    pmatrix.fs(bweim.list, t=c(5,10), trans=tmat)
   
    pmat1 <- pmatrix.fs(bweic.list, t=c(5,10), trans=tmat, newdata=list(x=-1))
    pmat3 <- pmatrix.fs(bweic.list, t=c(5,10), trans=tmat, newdata=list(x=c(-1,-1,-1)))
    expect_equal(pmat1[[1]], pmat3[[1]])

    expect_error(pmatrix.fs(bweic.list, t=c(5,10), trans=tmat, newdata=list(x=c(-1,-1))),
                 "must either have one row, or one row for each of the 3 allowed transitions")

    pmatrix.fs(bln.markov, t=5, trans=tmat)

    totlos.fs(bweim.list, t=5, trans=tmat)
    totlos.fs(bln.markov, t=5, trans=tmat)
})

test_that("list and non-list format give same estimates", { 
    expect_equal(
        as.numeric(pars.fmsm(bwei, trans=tmat)[[1]]),
        as.numeric(pars.fmsm(bwei.list, trans=tmat)[[1]]), tolerance=1e-04)
    bexpci <- flexsurvreg(Surv(years, status) ~ trans*x, data=bosms3, dist="exp")
    expect_equal(
        as.numeric(pars.fmsm(bexpci, newdata=list(x=1), trans=tmat)[[2]]),
        as.numeric(pars.fmsm(bexpc.list, newdata=list(x=1), trans=tmat)[[2]]), tolerance=1e-05)

    expect_equal(pmatrix.fs(bwei, trans=tmat), pmatrix.fs(bwei.list, trans=tmat), tolerance=1e-04)
    expect_equal(totlos.fs(bwei, trans=tmat), totlos.fs(bwei.list, trans=tmat), tolerance=1e-04)

    expect_equal(msfit.flexsurvreg(bwei, trans=tmat, t=1:10, variance=FALSE), 
                 msfit.flexsurvreg(bwei.list, trans=tmat, t=1:10, variance=FALSE), tolerance=1e-04)
    expect_equal(msfit.flexsurvreg(bexpci, newdata=list(x=1), trans=tmat, t=1:10, variance=FALSE),
                 msfit.flexsurvreg(bexpc.list, newdata=list(x=1), trans=tmat, t=1:10, variance=FALSE), tolerance=1e-05)
})



test_that("different model families for different transitions", { 

blist <- vector(3, mode="list")
blist[[1]] <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==1), data=bosms3, dist="exp")
blist[[2]] <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==2), data=bosms3, dist="gamma")
blist[[3]] <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==3), data=bosms3, dist="weibull")
expect_equal(pmatrix.fs(blist, t=c(5,10), trans=tmat)$`5`[1,2], 0.2423839, tolerance=1e-04)
expect_equal(totlos.fs(blist, trans=tmat)[1,2], 0.08277735, tolerance=1e-04)

blistx <- vector(3, mode="list")
blistx[[1]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==1), data=bosms3, dist="exp")
blistx[[2]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==2), data=bosms3, dist="gamma")
blistx[[3]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==3), data=bosms3, dist="weibull")
expect_equal(pmatrix.fs(blistx, t=c(5,10), trans=tmat, newdata=list(x=1))$`5`[1,2], 0.2291726, tolerance=1e-04)
expect_equal(totlos.fs(blistx, trans=tmat, newdata=list(x=1))[1,2],0.08698692 , tolerance=1e-04)

    
})
