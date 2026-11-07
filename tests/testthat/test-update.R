d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
fit <- jabba3(d, osa = FALSE)

test_that("update() without changes reproduces the fit", {
  f <- update(fit)
  expect_equal(f$par, fit$par, tolerance = 1e-6)
  expect_identical(f$settings$osa, FALSE)
  expect_identical(f$run, fit$run)
})

test_that("changing r keeps derived priors unless recompute = TRUE", {
  f1 <- update(fit, priors = list(r = c(0.3, 0.5)))
  expect_equal(f1$priors$r[["median"]], 0.3)
  expect_equal(f1$priors$sdb, fit$priors$sdb)
  expect_equal(f1$priors$kdyn, fit$priors$kdyn)
  f2 <- update(fit, priors = list(r = c(0.3, 0.5)), recompute = TRUE)
  expect_equal(f2$priors$sdb[["median"]], 0.80*0.5*0.3^0.45)
  expect_false(isTRUE(all.equal(f2$priors$kdyn, fit$priors$kdyn)))
})

test_that("dependent priors can be set explicitly", {
  f <- update(fit, priors = list(r = c(0.3, 0.5), sdb = c(0.15, 0.3)))
  expect_equal(unname(f$priors$sdb[1:2]), c(0.15, 0.3))
  expect_identical(f$priors$sdb.type, "lognormal")
})

test_that("model options and settings are updated", {
  f <- update(fit, model = list(lag = 3, rho = 0.75), run = "lag + AR(1)")
  expect_equal(f$model$lag, 3); expect_equal(f$model$rho, 0.75)
  expect_identical(f$run, "lag + AR(1)")
  f2 <- update(fit, model = list(shape = "schaefer"))
  expect_equal(f2$model$m, jb3_model(shape = "schaefer")$m)
})

test_that("objects passed whole are used as is, and bad names fail", {
  pr <- jb3_priors(d, r = c(0.25, 0.4))
  expect_equal(update(fit, priors = pr)$priors$r, pr$r)
  expect_error(update(fit, priors = list(rr = 1)), "unknown jb3_priors")
  expect_error(update(fit, model = list(lagg = 1)), "unknown jb3_model")
})

test_that("hand-edited models are patched directly with a warning", {
  f0 <- fit; f0$model$lag <- 2
  expect_warning(f <- update(f0, model = list(rho = 0.5)), "edited after jb3_model")
  expect_equal(f$model$lag, 2); expect_equal(f$model$rho, 0.5)
})

test_that("workflow fits rerun the workflow from the original inputs", {
  fw <- jabba3(d, workflow = "tier1", amat = 3, osa = FALSE)
  f <- update(fw, priors = list(r = c(0.3, 0.5)))
  expect_identical(f$settings$workflow, "tier1")
  expect_equal(f$model$lag, fw$model$lag)
  expect_length(f$steps, 2)
  expect_equal(f$steps[[1]]$model$lag, 0)
})
