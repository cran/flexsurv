if (!identical(Sys.getenv("NOT_CRAN"), "true")) return()

## simulate events following hospital 
n <- 1000
set.seed(1)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
events <- c("icu","death","discharge")
pbase <- c(0.2, 0.3, 0.5)
event <- numeric(n)
for (i in 1:n){
  p <- pmnlogit(qmnlogit(pbase) + 2*x[i] + 3*y[i])
  event[i] <- sample(events, size=1, prob=p, replace=TRUE)
}
t <- numeric(n)
t[event=="death"] <- rgamma(sum(event=="death"), 2.5, 1.2)
t[event=="discharge"] <- rgamma(sum(event=="discharge"), 3.5, 0.6)
t[event=="icu"] <- rgamma(sum(event=="icu"), 1, 3.2)
cens <- as.numeric(t > 3)
t[cens] <- 3
status <- 1 - cens
dat <- data.frame(t, status, x, y, event)

## model for event following hospital
fhosp <-   flexsurvmix(Surv(t, status) ~ x, pformula = ~x + y,
                       data=dat, event=event, 
                       dists=c("gamma","gamma","gamma"))

## simulate events following ICU
nicu <- sum(dat$event=="icu")
picu <- c(0.4, 0.6)
set.seed(1)
evicu <- sample(c("death","discharge"), size=nicu, prob=picu, replace=TRUE)
ti <- numeric(nicu)
ti[evicu=="death"] <- rgamma(sum(evicu=="death"), 1.5, 1)
ti[evicu=="discharge"] <- rgamma(sum(evicu=="discharge"), 0.5, 3)
censi <- as.numeric(ti > 1)
ti[censi] <- 1
statusi <- 1 - censi
dati <- data.frame(ti, statusi, evicu)

## model for event following ICU
ficu <- flexsurvmix(Surv(ti, statusi) ~ 1, data=dati, event=evicu, 
                       dists=c("gamma","gamma"))

## Construct multi-state model object
fm <- fmixmsm("hospital"=fhosp, "icu"=ficu)
nd <- data.frame(x=c(0,0.02), y=c(0,0.01))

test_that("ppath_fmixmsm",{
  probh <- probs_flexsurvmix(fhosp)
  probi <- probs_flexsurvmix(ficu)
  pp <- ppath_fmixmsm(fm)
  expect_equal(
    probh$val[probh$event=="icu"] * probi$val[probi$event=="death"],
    pp$val[pp$pathway=="hospital-icu-death"]
  )
  ppath_fmixmsm(fm, final=TRUE)
  probh <- probs_flexsurvmix(fhosp,newdata=nd)
  probi <- probs_flexsurvmix(ficu,newdata=nd)
  pp <- ppath_fmixmsm(fm, newdata=nd)
  expect_equal(
    probh$val[probh$event=="icu" & probh$x==0.02 & probh$y==0.01] * 
      probi$val[probi$event=="death" & probi$x==0.02 & probi$y==0.01],
    pp$val[pp$pathway=="hospital-icu-death" & pp$x==0.02 & pp$y==0.01]
  )
  ppath_fmixmsm(fm, newdata=nd, final=TRUE)
  expect_true(is.numeric(ppath_fmixmsm(fm, B=3)$lower))
  expect_true(is.numeric(ppath_fmixmsm(fm, final=TRUE, B=3)$lower))
  expect_true(is.numeric(ppath_fmixmsm(fm, newdata=nd, B=3)$lower))
  expect_true(is.numeric(ppath_fmixmsm(fm, newdata=nd, final=TRUE, B=3)$lower))
})

test_that("meanfinal_fmixmsm",{
  meanfinal_fmixmsm(fm)
  meanfinal_fmixmsm(fm, final=TRUE)
  meanfinal_fmixmsm(fm, newdata=nd)
  meanfinal_fmixmsm(fm, newdata=nd, final=TRUE)
  expect_true(is.numeric(meanfinal_fmixmsm(fm, B=3)$lower))
  expect_true(is.numeric(meanfinal_fmixmsm(fm, final=TRUE, B=3)$lower))
  expect_true(is.numeric(meanfinal_fmixmsm(fm, newdata=nd, B=3)$lower))
  expect_true(is.numeric(meanfinal_fmixmsm(fm, newdata=nd, final=TRUE, B=3)$lower))
})

if (covr::in_covr()){
  test_that("qfinal_fmixmsm",{
    expect_error({
      qfinal_fmixmsm(fm, newdata=nd)
      qfinal_fmixmsm(fm, newdata=nd, final=TRUE)
      qfinal_fmixmsm(fm, newdata=nd, probs=c(0.25, 0.75))
      qfinal_fmixmsm(fm, newdata=nd, probs=c(0.25, 0.75), final=TRUE)
      
      qfinal_fmixmsm(fm, newdata=nd, B=10)
      qfinal_fmixmsm(fm, newdata=nd, B=10, final=TRUE)
      qfinal_fmixmsm(fm, newdata=nd, n=100, B=10)
      qfinal_fmixmsm(fm, newdata=nd, n=100, B=10, final=TRUE)
    }, NA)
  })
}

