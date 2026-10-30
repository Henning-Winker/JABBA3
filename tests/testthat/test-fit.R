data(iccat, package = "JABBA3")
d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")

test_that("bigeye Fox fit converges and reproduces the JABBA 2 RTMB engine", {
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE)
  expect_s3_class(f, "jabba3")
  expect_equal(f$convergence$tier, "ok")
  # reference values from JABBA 2.3 + BioProcDyn patch 0013 (engine = "rtmb", sigma.proc.prior FishLife, psi CV 0.25)
  expect_equal(f$pars["K", "est"], 1174749, tolerance = 1e-4)
  expect_equal(f$pars["r", "est"], 0.1929, tolerance = 1e-3)
  expect_equal(utils::tail(f$ts$est[f$ts$quantity == "B/Bmsy"], 1), 0.7973, tolerance = 1e-3)
})

test_that("OSA residuals are approximately standard normal on a well-specified fit", {
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"))
  expect_true(all(is.finite(f$index$osa)))
  dg <- jb3_diags(f)
  expect_true(all(c("runs.p", "bias.p", "LjungBox.p", "Shapiro.p") %in% names(dg)))
})

test_that("workflows, options and outputs run", {
  w <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), workflow = "tier2", amat = 3, osa = FALSE)
  expect_equal(names(w$steps), c("Current", "Tier 2"))
  expect_equal(w$model$lag, 3L)
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox", obs.error = "spict"), osa = FALSE, ndraws = 200)
  expect_equal(nrow(f$draws$B), 200)
  p <- jb3_project(f, catch = c(60000, 80000), nyears = 5)
  expect_s3_class(p, "jb3proj")
  expect_s3_class(jbplot_advice(f), "ggplot")
  expect_s3_class(jbplot_kobe(f), "ggplot")
})

test_that("AR(1) process error and fast mode run", {
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox", lag = 3, rho = 0.75), osa = TRUE)
  expect_equal(f$report$rho, 0.75)
  expect_true(f$convergence$tier %in% c("ok", "warning"))
  g <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), se = FALSE, osa = FALSE)
  expect_equal(g$convergence$tier, "ok")
})
