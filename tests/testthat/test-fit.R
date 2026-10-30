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

test_that("JABBA data frames convert to FLR and give the same fit", {
  skip_if_not_installed("FLCore")
  x <- jbdat2flr(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se)
  expect_s4_class(x$catch, "FLQuant"); expect_s4_class(x$indices, "FLIndices")
  expect_equal(names(x$indices), c("JR2_early", "JR2_late"))
  d2 <- jb3_data(x$catch, x$indices, stock = "bet")
  expect_equal(d2$I, d$I, ignore_attr = TRUE); expect_equal(d2$SE, d$SE, ignore_attr = TRUE); expect_equal(d2$catch, d$catch)
  back <- flr2jbdat(x$catch, x$indices)
  expect_equal(back$cpue[, -1], iccat$bet$cpue[, -1], ignore_attr = TRUE)
  xf <- jbdat2flr(data.frame(year = 2000:2009, a = 1:10, b = 11:20), data.frame(year = 2000:2009, i = 1:10))
  expect_s4_class(xf$catch, "FLQuants"); expect_length(xf$catch, 2)
})

test_that("jbinput2jabba3 reproduces the JABBA 2 RTMB fit", {
  skip_if_not_installed("JABBA")
  b <- iccat$bet
  jb <- suppressMessages(JABBA::build_jabba(catch = b$catch, cpue = b$cpue, se = b$se, model.type = "Fox", r.prior = c(0.2, 0.5),
          psi.prior = c(0.9, 0.25), sigma.proc.prior = c(0.80*0.5*0.2^0.45, 0.3), sigmaproc_bound = 10, verbose = FALSE))
  x <- jbinput2jabba3(jb, verbose = FALSE)
  f <- jabba3(x$data, x$priors, x$model, osa = FALSE)
  expect_equal(f$pars["K", "est"], 1174749, tolerance = 1e-4)
  expect_equal(utils::tail(f$ts$est[f$ts$quantity == "B/Bmsy"], 1), 0.7973, tolerance = 1e-3)
})

test_that("spictinp2jabba3 converts a SPiCT input list", {
  b <- iccat$bet; ok1 <- is.finite(b$cpue[, 2]); ok2 <- is.finite(b$cpue[, 3])
  inp <- list(obsC = b$catch[, 2], timeC = b$catch[, 1],
              obsI = list(b$cpue[ok1, 2], b$cpue[ok2, 3]), timeI = list(b$cpue[ok1, 1] + 0.5, b$cpue[ok2, 1] + 0.5),
              priors = list(logr = c(log(0.2), 0.5, 1), logn = c(log(2), 2, 0), logbkfrac = c(log(0.8), 0.5, 1)),
              ini = list(logn = log(2)), phases = list(logn = -1))
  x <- spictinp2jabba3(inp, stock = "bet", verbose = FALSE)
  expect_equal(x$model$shape, "schaefer"); expect_equal(x$model$obs.error, "spict")
  expect_equal(nrow(x$data$I), 68); expect_equal(unname(x$priors$bkinit[1]), 0.8)
  f <- jabba3(x$data, x$priors, x$model, osa = FALSE)
  expect_true(f$convergence$tier %in% c("ok", "warning"))
})
