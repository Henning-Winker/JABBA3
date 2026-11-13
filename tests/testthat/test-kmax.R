test_that("K.max (deprecated): r-based or multiple-of-max-catch cap still works, binds only above the cap", {
  data(iccat); bet <- iccat$bet; mc <- max(bet$catch[[2]])
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  expect_equal(unname(jb3_priors(d, K = "jabba")$K.max[3]), 0)                  # two-sided JABBA prior: no cap
  expect_message(pk <- jb3_priors(d, r = c(0.2, 0.5), K = "jabba", K.max = "r"), "deprecated")
  expect_equal(unname(pk$K.max[1]), 7.6*0.2^-0.587, tolerance = 1e-6)
  f0 <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5), K = "jabba"), jb3_model(shape = "fox"), osa = FALSE)
  f1 <- jabba3(d, pk, jb3_model(shape = "fox"), osa = FALSE)
  expect_equal(f0$pars["K", "est"], f1$pars["K", "est"], tolerance = 1e-3)     # K well below the cap: no effect
  f2 <- suppressMessages(jabba3(d, jb3_priors(d, r = c(0.2, 0.5), K.max = c(4, 0.05)), jb3_model(shape = "fox"), osa = FALSE))
  expect_lt(f2$pars["K", "est"]/mc, 4.5)                                         # tight cap below the estimate binds
})

test_that("K prior shortcuts: half8 default, jabba (8 x max catch, CV 1), jabbaN, smoothed max catch", {
  data(iccat); bet <- iccat$bet; mc <- max(bet$catch[[2]])
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  p <- jb3_priors(d); s1 <- sqrt(log(2))
  expect_equal(unname(p$K[3]), 0); expect_equal(unname(p$K.max[1:3]), c(8*exp(-s1^2/2), s1, 1))
  expect_equal(unname(jb3_priors(d, K = "jabba")$K[1:3]), c(8*mc, s1, 1), tolerance = 1e-8)
  expect_equal(unname(jb3_priors(d, K = "jabba4")$K[1]), 4*mc, tolerance = 1e-8)
  ps <- jb3_priors(d, K = "jabba4", maxC = "smooth")
  expect_lt(ps$maxC.value, mc); expect_equal(unname(ps$K[1]), 4*ps$maxC.value, tolerance = 1e-8)
  expect_error(jb3_priors(d, K = "foo"), "jabba")
})

test_that("K plausibility warning and K-profile starts in the best-of-starts search", {
  data(rjc)                                                   # thornback ray 4a: K ~ 50 x max catch, above the r-based bound
  x <- suppressMessages(spictinp2jabba3(rjc$inp, stock = "rjc"))
  expect_equal(unname(x$priors$K.max[3]), 0)
  expect_warning(f <- jabba3(x$data, x$priors, x$model, osa = FALSE), "plausibility bound")
  expect_gt(f$pars["K", "est"]/max(x$data$catch, na.rm = TRUE), 40)  # not pulled to the bound
})
