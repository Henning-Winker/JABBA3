test_that("tables and F-based projections", {
  data(iccat); bet <- iccat$bet
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE)
  pt <- jb3_partable(f); expect_true(all(c("K", "r", "Fmsy", "MSY", "rho") %in% pt$par))
  expect_equal(ncol(jb3_partable(list(a = f, b = f), format = "wide")), 3)
  ts <- jb3_tstable(f); expect_equal(nrow(ts), length(f$years)); expect_true("F/Fmsy.uci" %in% names(ts))
  pj <- jb3_project(f, nyears = 5, ndraws = 200)
  expect_equal(pj$scenarios$type[1], "F/Fmsy")
  tb <- jb3_projtable(pj); expect_equal(nrow(tb), 5)
  expect_equal(tb[["P(F<=Fmsy) 2018"]][tb$scenario == "F/Fmsy = 1"], 1)
  km <- jb3_kobematrix(pj); expect_true(all(km[, -1] >= 0 & km[, -1] <= 1))
  pc <- jb3_project(f, catch = 1e7, nyears = 3, ndraws = 100)            # catch capped: F stays finite
  expect_true(all(is.finite(pc$iters[[1]]$`F/Fmsy`)))
})
