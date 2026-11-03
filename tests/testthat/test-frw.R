test_that("missing catch is estimated with an F random walk (catch.gaps = 'estimate')", {
  skip_if_not_installed("RTMB")
  data(iccat); bet <- iccat$bet
  expect_error(jb3_data(bet$catch[!bet$catch$year %in% 1990:1994, ], bet$cpue, bet$se), "catch.gaps")
  d <- suppressMessages(jb3_data(bet$catch[!bet$catch$year %in% 1990:1994, ], bet$cpue, bet$se, catch.gaps = "estimate"))
  expect_equal(sum(!d$catch.obs), 5)
  f <- suppressWarnings(jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE))
  expect_true(f$convergence$tier %in% c("ok", "warning"))
  cc <- f$ts[f$ts$quantity == "Catch", ]; tru <- bet$catch[[2]][bet$catch$year %in% 1990:1994]
  ratio <- cc$est[cc$year %in% 1990:1994]/tru
  expect_true(all(ratio > 0.7 & ratio < 1.4))
  expect_true("sdf" %in% rownames(f$pars))
  # catch starting after the index: model starts at the first index year
  d2 <- suppressMessages(jb3_data(bet$catch[bet$catch$year >= 1970, ], bet$cpue, bet$se, catch.gaps = "estimate"))
  expect_equal(min(d2$years), 1959)
  # complete catch: identical to the default path
  d0 <- jb3_data(bet$catch, bet$cpue, bet$se); d1 <- jb3_data(bet$catch, bet$cpue, bet$se, catch.gaps = "estimate")
  f0 <- suppressWarnings(jabba3(d0, jb3_priors(d0, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE))
  f1 <- suppressWarnings(jabba3(d1, jb3_priors(d1, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE))
  expect_equal(f0$opt$objective, f1$opt$objective)
})
