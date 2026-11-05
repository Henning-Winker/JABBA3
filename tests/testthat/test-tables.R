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
  pk <- jb3_project(f, catch = c(5e4, 7e4, 9e4), nyears = 4, ndraws = 200, imp.years = 1)
  k2 <- jb3_kobematrix(pk, what = "all", percent = TRUE)
  expect_s3_class(k2, "jb3k2sm"); expect_equal(names(k2[[1]])[1], "TAC"); expect_equal(nrow(k2$`B>Bmsy`), 3)
  expect_equal(length(unique(k2$`B>Bmsy`[["2018"]])), 1)                  # common random numbers: same intermediate year
  expect_true(all(diff(k2$`F<Fmsy`[["2021"]]) <= 0))                       # higher TAC, lower P(F < Fmsy)
  expect_true(all(k2$`F<Fmsy & B>Bmsy`[, -1] <= k2$`B>Bmsy`[, -1]))
  kc <- jb3_kobematrix(list(pk, pk), years = 2019:2020); expect_equal(kc, jb3_kobematrix(pk, years = 2019:2020))
  expect_error(jb3_kobematrix(pk, years = 2030))
  ki <- jb3_kobematrix_iotc(pk, years = c(2019, 2021), ref.catch = 7e4, blim = 0.4, flim = 1.4)
  expect_s3_class(ki, "jb3k2iotc"); expect_equal(nrow(ki), 8); expect_equal(names(ki)[4], "71% (50,000t)")
  expect_true(all(ki[ki$refpt == "B < Blim", -(1:3)] <= ki[ki$refpt == "B < Bmsy", -(1:3)]))
  expect_s3_class(jbplot_kobematrix(ki), "ggplot"); expect_s3_class(jbplot_kobematrix(ki, colours = "kobe"), "ggplot")
  expect_output(print(ki), "limit reference points")
  kb <- jb3_kobematrix(pk, what = "all", percent = TRUE, bref = c(Btrigger = 0.5))
  expect_equal(names(kb), c("F<Fmsy", "B>Btrigger", "F<Fmsy & B>Btrigger"))
  kl <- jb3_kobematrix(pk, percent = TRUE, bref = c(Blim = 0.3), fref = NULL, risk = TRUE)
  expect_equal(attr(kl, "condition"), "B<Blim")
  expect_equal(kl[, -1], 100 - jb3_kobematrix(pk, percent = TRUE, bref = c(Blim = 0.3), fref = NULL)[, -1])
  expect_error(jb3_kobematrix(pk, what = "F", fref = NULL))
  expect_s3_class(jbplot_kobematrix(pk, bref = c(Blim = 0.3), fref = NULL, risk = TRUE), "ggplot")
  pr <- jb3_project(f, scenarios = list(Fmsy = c(Fmsy = 1), Fcur = c(Fsq = 1), "0.8Fmsy" = c(Fmsy = 0.8), TAC = c(Catch = 6e4)),
                    nyears = 3, ndraws = 100)
  kr <- jb3_kobematrix(pr); expect_equal(kr$scenario, c("Fmsy", "Fcur", "0.8Fmsy", "TAC"))
  expect_equal(jb3_project(f, fmsy = c(a = 1, b = 0.5), nyears = 2, ndraws = 50)$scenarios$label, c("a", "b"))
  expect_s3_class(jbplot_kobematrix(pr), "ggplot")
  expect_s3_class(jbplot_kobematrix(pk), "ggplot"); expect_s3_class(jbplot_kobematrix(k2, what = "green", colours = "kobe"), "ggplot")
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

test_that("jabbaplot3 panels", {
  skip_if_not_installed("patchwork")
  data(iccat); bet <- iccat$bet
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"))
  expect_s3_class(plot(f), "patchwork")
  expect_s3_class(jabbaplot3(list(a = f, b = f), panels = c("B/Bmsy", "kobe", "production", "procdev"), ncol = 2), "patchwork")
  expect_error(jabbaplot3(f, panels = "nope"), "unknown panels")
})

test_that("rho guidance", {
  expect_equal(jb3_rho_gt(8), round(exp(-2.2/8), 2))
  expect_equal(nrow(jb3_rho_table()), 4)
})

test_that("Pella shape with Bmsy/K below 1/e", {
  m <- jb3_model(BmsyK = 0.36)$m
  expect_lt(m, 1); expect_equal(m^(-1/(m - 1)), 0.36, tolerance = 1e-6)
})
