test_that("output columns follow FLR naming: year, quant, season, iter", {
  data(iccat); b <- iccat$bum
  d <- jb3_data(b$catch, b$cpue, b$se)
  f <- jabba3(d, jb3_priors(d, r = c(0.1, 0.3)), jb3_model(), osa = FALSE)
  expect_true(all(c("year", "quant", "est", "lci", "uci") %in% names(f$ts)))
  expect_false(any(c("quantity", "Year") %in% names(f$ts)))
  p <- jb3_project(f, catch = 2000, nyears = 3, ndraws = 50)
  expect_true("quant" %in% names(p$ts))
  expect_true(!is.null(p$iter) && is.null(p$iters))
  r <- jb3_retro(f, n = 1)
  expect_true("quant" %in% names(r$rho))
})
