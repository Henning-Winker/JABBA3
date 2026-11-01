test_that("jb3_contraction is zero against itself and positive with an index", {
  skip_on_cran()
  data(iccat); bet <- iccat$bet
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  pr <- jb3_priors(d, r = c(0.2, 0.5))
  f <- jabba3(d, pr, jb3_model(shape = "fox"), osa = FALSE, ndraws = 500)
  ct0 <- jb3_contraction(f, prior = f, quantities = c("K", "r", "B/Bmsy"))
  expect_true(all(abs(ct0$contraction) < 1e-8))
  expect_s3_class(ct0, "jb3contraction")
  expect_s3_class(jbplot_contraction(ct0), "ggplot")
})
