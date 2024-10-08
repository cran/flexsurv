test_that("check tidy", {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")
  expect_equal(
    tidy(fitw, transform = 'baseline.real')$estimate,
    coef(fitw)
  )
})

test_that("check tidy against survreg", {
  fit_fs <- flexsurvreg(Surv(futime, fustat) ~ ecog.ps,
                      data = ovarian, dist = "weibull")
  fit_sr <- survreg(Surv(futime, fustat) ~ ecog.ps,
                    data = ovarian, dist = "weibull")
  expect_equal(
    object = as.data.frame(tidy(fit_fs)[3, ]),
    expected = as.data.frame(broom::tidy(fit_sr)[2, ]),
    tolerance = 1e-6
  )
})

test_that("check glance", {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")
  gl <- glance(fitw)
  expect_equal(gl$N, fitw$N)
  expect_equal(gl$events, fitw$events)
  expect_equal(gl$trisk, fitw$trisk)
  expect_equal(gl$df, fitw$npars)
  expect_equal(gl$logLik, fitw$loglik)
  expect_equal(gl$AIC, fitw$AIC)
  expect_equal(gl$BIC, BIC(fitw))
})

test_that("check augment", {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                      dist = "weibull")
  expect_equal(augment(fitw)$.pred_time, predict(fitw)$.pred_time)
  expect_equal(augment(fitw)$.resid, residuals(fitw))
})
