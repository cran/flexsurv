fit <- flexsurvreg(formula = Surv(futime, fustat) ~ rx, data = ovarian, dist="weibull")
nd <- data.frame(rx=1:2)

test_that("simulate.flexsurvreg",{
  sim <- simulate(fit)
  expect_true(nrow(sim)==nrow(ovarian))
  sim <- simulate(fit, newdata=nd)
  expect_true(all(sim$time_1 > 0))
  sim <- simulate(fit, seed=1002, newdata=nd)
  expect_equal(sim$time_1, c(575, 2392), tolerance=1)
  sim <- simulate(fit, seed=1002, newdata=nd, nsim=5)
  expect_true(all(sim$time_2 > 0))
  sim <- simulate(fit, seed=1002, newdata=nd, nsim=5, censtime = 1000, tidy=TRUE)
  expect_equal(sim$event[sim$time==1000], rep(0, 5))  
  sim <- simulate(fit, seed=1002, newdata=nd, nsim=5, censtime = c(500,1000), tidy=TRUE)
  expect_equal(sim$event, c(0,0,0,1,1,0,0,0,1,0))
  sim <- simulate(fit, seed=1002, newdata=nd, start=500, tidy=TRUE)
  expect_true(all(sim$time > 500))
  sim <- simulate(fit, seed=1002, newdata=nd, start=c(500,700), tidy=TRUE)
  expect_true(all(sim$time > 500))
})

test_that("simulate.flexsurvreg with left truncation",{
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ rx, data = ovarian, dist="weibull")
  nd <- ovarian
  sim <- simulate(fit, seed=1003, newdata=nd, nsim = 20, start = nd$futime)
  expect_true(all(sim[,1:20] > nd$futime))
})
