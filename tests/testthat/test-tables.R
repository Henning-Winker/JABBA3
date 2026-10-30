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

test_that("FLStockR and procerr options", {
  skip_if_not_installed("FLCore")
  data(iccat); bet <- iccat$bet
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox", rho = 0.75), osa = FALSE)
  s <- jb3_FLStockR(f, ndraws = 50)
  expect_s4_class(s, "FLStockR"); expect_equal(dim(s@stock)[6], 50)
  expect_equal(unname(c(s@refpts["Blim"]/s@refpts["Bmsy"])), 0.3)
  for(pe in c("ar1", "mean", "none")) expect_s3_class(jb3_project(f, fmsy = 1, nyears = 3, ndraws = 50, procerr = pe), "jb3proj")
})
