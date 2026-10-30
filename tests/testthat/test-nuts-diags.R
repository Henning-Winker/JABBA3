test_that("NUTS diagnostics and plots", {
  skip_if_not_installed("tmbstan"); skip_on_cran()
  data(iccat); bet <- iccat$bet
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  fn <- suppressWarnings(jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), engine = "nuts",
                                nuts = list(chains = 2, iter = 300, warmup = 150)))
  dg <- jb3_nuts_diags(fn)
  expect_true(all(c("Rhat", "ess_bulk", "ess_tail") %in% names(dg$pars)))
  expect_equal(nrow(dg$sampler), 2)
  for(t in c("trace", "rank", "dens", "acf", "pairs", "rhat", "energy")) expect_s3_class(jbplot_nuts(fn, t, pars = c("logK", "logr", "logq")[if(t %in% c("rhat", "energy")) 0 else 1:3]), "ggplot")
  expect_error(jb3_nuts_diags(jabba3(d, jb3_priors(d), osa = FALSE)), "engine")
})
