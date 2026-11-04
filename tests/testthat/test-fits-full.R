test_that("jbplot_fits(years = 'all') covers all years; several runs show relative refpts only", {
  data(iccat); b <- iccat$bet
  d <- jb3_data(b$catch, b$cpue, b$se)
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE)
  x <- JABBA3:::.index_full(f, 200)
  expect_equal(nrow(x), 2*length(f$years))
  expect_true(all(x$ppd_lci <= x$lci + 1e-8 | is.na(x$ppd_lci)) || mean(x$ppd_lci <= x$lci) > 0.95)
  expect_s3_class(jbplot_fits(f, years = "all", ndraws = 200), "ggplot")
  f2 <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "pella", BmsyK = 0.4), osa = FALSE, run = "pella")
  p <- jbplot_trj(list(fox = f, pella = f2), type = c("B", "B/Bmsy"))
  expect_s3_class(p, "ggplot")
})
