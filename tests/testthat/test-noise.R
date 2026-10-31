test_that("jb3_splitindex splits with and without overlap", {
  I <- data.frame(year = 2001:2010, a = 1:10, b = 10:1); se <- data.frame(year = 2001:2010, a = 0.2, b = 0.3)
  x <- jb3_splitindex(I, se, "a", at = 2006)
  expect_equal(names(x$I), c("year", "a.old", "a.new", "b"))
  expect_true(all(is.na(x$I$a.old[6:10])) && all(is.na(x$I$a.new[1:5])))
  x1 <- jb3_splitindex(I, se, "a", at = 2006, overlap = 2006:2007, new = c(60, 70))
  expect_equal(x1$I$a.old[6:7], 6:7); expect_equal(x1$I$a.new[6:7], c(60, 70))
  expect_equal(x1$se$a.new[1:5], rep(NA_real_, 5))
})
test_that("jb3_noise returns observation and loess CVs", {
  set.seed(1); y <- 1991:2010
  I <- data.frame(year = y, a = exp(cumsum(rnorm(20, 0, 0.1)) + rnorm(20, 0, 0.2)))
  z <- jb3_noise(I, data.frame(year = y, a = 0.2))
  expect_equal(z$obs.cv, 0.2); expect_true(z$loess.cv > 0.05 && z$loess.cv < 0.5)
})
