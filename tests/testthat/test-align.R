test_that("inputs of different lengths are aligned to the catch years", {
  data(iccat); bet <- iccat$bet
  d0 <- jb3_data(bet$catch, bet$cpue, bet$se)
  # catch shorter than the indices: index years before the catch are not used
  expect_message(d1 <- jb3_data(bet$catch[bet$catch$year >= 1970, ], bet$cpue, bet$se), "not used")
  expect_equal(d1$years, 1970:2017)
  expect_equal(d1$I, d0$I[as.character(1970:2017), ])
  # rows in any order, se for some years only
  set.seed(1)
  expect_equal(jb3_data(bet$catch[sample(nrow(bet$catch)), ], bet$cpue[sample(nrow(bet$cpue)), ], bet$se)$I, d0$I)
  d2 <- jb3_data(bet$catch, bet$cpue, bet$se[bet$se$year >= 2000, ])
  expect_true(all(is.na(d2$SE[as.character(1950:1999), ])))
  expect_equal(d2$SE[as.character(2000:2017), ], d0$SE[as.character(2000:2017), ])
  # se matched by name (reordered, subset)
  se <- bet$se; names(se) <- names(bet$cpue)
  expect_equal(jb3_data(bet$catch, bet$cpue, se[, c(1, 3, 2)])$SE, jb3_data(bet$catch, bet$cpue, se)$SE)
  d3 <- jb3_data(bet$catch, bet$cpue, se[, c(1, 3)])
  expect_true(all(is.na(d3$SE[, 1]))); expect_equal(d3$SE[, 2], d0$SE[, 2])
  expect_error(jb3_data(bet$catch, bet$cpue, data.frame(year = se$year, JR2_late = se[[3]], other = 1)), "named as the indices")
  # catch padded with NA years is trimmed; gaps inside are an error
  cp <- rbind(data.frame(year = 1945:1949, catch = NA), setNames(bet$catch[, 1:2], c("year", "catch")))
  expect_equal(jb3_data(cp, bet$cpue, bet$se)$years, d0$years)
  expect_error(jb3_data(bet$catch[bet$catch$year != 1990, ], bet$cpue, bet$se), "catch missing in years 1990")
  # fleets: NA counts as no catch unless all fleets are missing
  cf <- data.frame(year = bet$catch$year, f1 = bet$catch[[2]]*0.6, f2 = ifelse(bet$catch$year < 1970, NA, bet$catch[[2]]*0.4))
  d4 <- jb3_data(cf, bet$cpue, bet$se)
  expect_equal(d4$catch[d4$years >= 1970], d0$catch[d0$years >= 1970])
  expect_equal(nrow(jb3_align(d0)), 3)
})

test_that("FLR inputs of different lengths are aligned", {
  skip_if_not_installed("FLCore")
  suppressMessages(library(FLCore))
  data(iccat); bet <- iccat$bet
  x <- jbdat2flr(bet$catch, bet$cpue, bet$se)
  i1 <- window(x$indices[[1]], start = 1960, end = 1990); i2 <- window(x$indices[[2]], start = 1975, end = 2017)
  d <- jb3_data(x$catch, FLIndices(a = i1, b = i2))
  expect_equal(d$years, 1950:2017)
  expect_message(d2 <- jb3_data(window(x$catch, start = 1965), FLIndices(a = i1, b = i2)), "not used")
  a <- jb3_align(d2); expect_equal(a$n_used[a$series == "a"], a$n[a$series == "a"] - 5)
})
