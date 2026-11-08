test_that("jbplot_advice() takes FLStock(R), FLStocks and mixed lists", {
  skip_if_not_installed("FLCore")
  d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
  fit <- jabba3(d, osa = FALSE, run = "JABBA3")
  stk <- jb3_FLStockR(fit, ndraws = 200)
  st1 <- jb3_FLStockR(fit)                                     # single iteration: trajectory only
  f <- JABBA3:::.flr2fit(stk, "x")
  expect_equal(sort(unique(f$ts$quantity)), sort(c("B/Bmsy", "F/Fmsy", "Catch", "Catch/MSY")))
  expect_equal(f$ts$est[f$ts$quantity == "B/Bmsy"], fit$ts$est[fit$ts$quantity == "B/Bmsy"], tolerance = 0.1)
  expect_true(all(is.na(JABBA3:::.flr2fit(st1, "y")$ts$lci)))
  cl <- c(JABBA3 = "dodgerblue", Other = "black")
  for(p in list(jbplot_advice(stk), jbplot_advice(FLCore::FLStocks(JABBA3 = stk, Other = st1), colours = cl),
                jbplot_advice(list(JABBA3 = fit, Other = st1), colours = cl),
                jbplot_advice(jb3_FLStockR(fit, rel = FALSE), rel = FALSE)))
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_error(jbplot_advice(list(1, 2)), "FLStock")
})
