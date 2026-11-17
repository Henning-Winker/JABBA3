test_that("jbinput2jabba3() reconstructs a JABBA fit: data, SE2, priors, igamma proxy, JABBA estimates", {
  jf <- readRDS(test_path("fixtures", "jbfit_bet_diffuse_min.rds"))     # JABBA 2 Fox fit of ICCAT bigeye, igamma(0.001, 0.001)
  x <- suppressMessages(jbinput2jabba3(jf, verbose = FALSE))
  s <- jf$settings
  expect_equal(unname(x$data$catch), jf$inputseries$catch[[2]])
  se2 <- sqrt(as.matrix(s$SE2)); ok <- !is.na(as.matrix(jf$inputseries$cpue[, -1]))
  expect_equal(unname(x$data$SE[ok]), unname(se2[ok]))           # JABBA total input variance carried over
  expect_equal(x$model$se.min, 0)
  expect_equal(unname(x$priors$r[1:2]), as.numeric(s$r.pr[1:2]))
  expect_equal(unname(x$priors$K[1:3]), c(as.numeric(s$K.pr[1:2]), 1)); expect_equal(unname(x$priors$K.max[3]), 0)
  ps <- sqrt(jf$pars_posterior$sigma2)
  expect_equal(unname(x$priors$sdb[1:2]), c(median(ps), sd(log(ps))))              # diffuse igamma -> posterior proxy
  expect_equal(x$model$shape, "fox"); expect_true(all(c("B/Bmsy", "F/Fmsy", "B/K") %in% x$jabba$ts$quant))
  xf <- suppressMessages(jbinput2jabba3(jf, igamma.proxy = "fixed", verbose = FALSE)); expect_equal(xf$model$sdb.fixed, median(ps))
  xi <- suppressMessages(jbinput2jabba3(jf, igamma.proxy = "igamma", verbose = FALSE))
  expect_equal(xi$priors$sdb.type, "igamma"); expect_equal(xi$priors$igamma, c(0.001, 0.001))
  jf4 <- jf; jf4$settings$igamma <- c(4, 0.01)                                       # informative igamma: exact
  x4 <- suppressMessages(jbinput2jabba3(jf4, verbose = FALSE)); expect_equal(x4$priors$igamma, c(4, 0.01))
  f <- suppressWarnings(jabba3(x$data, x$priors, jb3_model(shape = "fox", se.min = 0, sdi.est = FALSE), osa = FALSE))
  b <- merge(f$ts[f$ts$quant == "B/Bmsy", c("year", "est")], x$jabba$ts[x$jabba$ts$quant == "B/Bmsy", c("year", "est")], by = "year")
  expect_lt(max(abs(log(b$est.x/b$est.y))), 0.1)                                    # JABBA3 mode within 10% of the JABBA medians
})

test_that("igamma(a, b) prior on the process variance: informative prior gives the mode near sqrt(b/a) without data", {
  data(iccat); bet <- iccat$bet; d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  p <- jb3_priors(d, sdb = "igamma", igamma = c(4, 0.01)); expect_equal(p$igamma, c(4, 0.01))
  expect_null(jb3_priors(d)$igamma)
  f <- suppressWarnings(jabba3(d, p, jb3_model(shape = "fox"), osa = FALSE))
  expect_gt(f$pars["sdb", "est"], 0.02); expect_lt(f$pars["sdb", "est"], 0.1)
})
