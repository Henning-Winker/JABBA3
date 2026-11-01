test_that("sdi.est = FALSE fixes the additional observation SD (JABBA sigma.est = FALSE)", {
  data(iccat); bet <- iccat$bet
  d <- jb3_data(bet$catch, bet$cpue, stock = "bet")
  f <- jabba3(d, jb3_priors(d), jb3_model(shape = "fox", se.min = 0.2, sdi.est = FALSE), osa = FALSE)
  expect_true(all(f$report$sdi < 1e-3))
})
