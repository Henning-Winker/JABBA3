test_that("K.max soft upper prior: r-based cap by default, binds only above the cap", {
  data(iccat); bet <- iccat$bet; mc <- max(bet$catch[[2]])
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  expect_equal(unname(jb3_priors(d)$K.max[3]), 1)                               # default "r"
  expect_equal(unname(jb3_priors(d, K.max = NULL)$K.max[3]), 0)
  pk <- jb3_priors(d, r = c(0.2, 0.5), K.max = "r")
  expect_equal(unname(pk$K.max[1]), 7.6*0.2^-0.587, tolerance = 1e-6)
  f0 <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5), K.max = NULL), jb3_model(shape = "fox"), osa = FALSE)
  f1 <- jabba3(d, pk, jb3_model(shape = "fox"), osa = FALSE)
  expect_equal(f0$pars["K", "est"], f1$pars["K", "est"], tolerance = 1e-3)     # K well below the cap: no effect
  f2 <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5), K.max = c(4, 0.05)), jb3_model(shape = "fox"), osa = FALSE)
  expect_lt(f2$pars["K", "est"]/mc, 4.5)                                         # tight cap below the estimate binds
})

test_that("K prior shortcuts: jabba (8 x max catch, CV 1), jabbaN, smoothed max catch", {
  data(iccat); bet <- iccat$bet; mc <- max(bet$catch[[2]])
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  expect_equal(unname(jb3_priors(d)$K[1:2]), c(8*mc, sqrt(log(2))), tolerance = 1e-8)
  expect_equal(unname(jb3_priors(d, K = "jabba")$K[1]), 8*mc, tolerance = 1e-8)
  expect_equal(unname(jb3_priors(d, K = "jabba4")$K[1]), 4*mc, tolerance = 1e-8)
  ps <- jb3_priors(d, K = "jabba4", maxC = "smooth")
  expect_lt(ps$maxC.value, mc); expect_equal(unname(ps$K[1]), 4*ps$maxC.value, tolerance = 1e-8)
  expect_error(jb3_priors(d, K = "foo"), "jabba")
})
