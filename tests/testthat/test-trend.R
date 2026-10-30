test_that("relative joint index and diagnostics", {
  data(iccat); b <- iccat$bum
  tr <- jb3_trend(b$cpue, b$se, fixed.obsE = 0.1)
  expect_s3_class(tr, "jb3trend"); expect_equal(tr$convergence$tier, "ok")
  expect_equal(tr$q[1], 1)
  expect_true(all(is.finite(tr$index$osa)))
  ji <- jb3_jointindex(tr)
  expect_s3_class(ji, "jb3jointindex"); expect_equal(nrow(ji), length(tr$years)); expect_equal(mean(ji$trend), 1)
  expect_true(all(is.finite(jb3_jointindex(tr, fill = TRUE)$index)))
  for(p in list(jbplot_fits(tr), jbplot_runs(tr, type = "ordinary"), jbplot_residuals(tr), jbplot_ppc(tr, ndraws = 200),
                jbplot_trend(tr), jbplot_jointindex(tr))) expect_s3_class(p, "ggplot")
  d <- jb3_data(b$catch, ji, stock = "bum")
  expect_equal(ncol(d$I), 1)
  f <- jabba3(d, jb3_priors(d, r = c(0.1, 0.2)), jb3_model(BmsyK = 0.36), osa = FALSE)
  expect_equal(f$convergence$tier, "ok")
  if(requireNamespace("FLCore", quietly = TRUE)){
    fi <- jb3_jointindex(tr, "FLIndexBiomass"); expect_s4_class(fi, "FLIndexBiomass")
    d2 <- jb3_data(b$catch, FLCore::FLIndices(joint = fi)); expect_equal(unname(d2$SE[, 1]), unname(d$SE[, 1]))
  }
})

test_that("census and catch reconstruction", {
  set.seed(3); yrs <- 1990:2020; tru <- sapply(c(100, 50, 30), function(a) a*exp(cumsum(c(0, rnorm(30, -0.01, 0.1)))))
  obs <- tru*exp(matrix(rnorm(length(tru), 0, 0.2), ncol = 3)); obs[sample(length(obs), 20)] <- NA
  tc <- jb3_trend(data.frame(year = yrs, A = obs[, 1], B = obs[, 2], C = obs[, 3]),
                  data.frame(year = yrs, A = 0.2, B = 0.2, C = 0.2), type = "census", variance.weighting = "fixed")
  expect_equal(nrow(tc$total), length(yrs)); expect_true(all(tc$total$est > 0))
  data(iccat); cc <- iccat$bet$catch; cc[c(5:8, 30:31), 2] <- NA
  cf <- jb3_catchfill(cc)
  expect_true(all(is.finite(cf$catch))); expect_equal(sum(cf$filled), 6)
  expect_equal(cf$catch[!cf$filled], cc[[2]][!cf$filled])
  expect_true(all(cf$se.log[cf$filled] > cf$se.log[!cf$filled][1]))
  d <- jb3_data(cf, iccat$bet$cpue, iccat$bet$se)
  expect_equal(length(d$catch.se), length(d$years))
  f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE)
  expect_equal(f$convergence$tier, "ok")
  expect_equal(length(f$setup$dat$sdc), length(d$years))
})

test_that("time block and multi-series catch reconstruction (Dentex)", {
  data(jara_examples); L <- jara_examples$dentex; L <- L[L$year >= 1976, ]
  cf <- jb3_catchfill(L, se = 0.001, timeblock = 1992)
  fit <- attr(cf, "fit")
  expect_equal(dim(fit$r), c(2L, 2L)); expect_true(all(fit$r[1, ] < 0) && all(fit$r[2, ] > 0))
  tab <- jb3_fill_table(cf)
  expect_equal(names(tab), c("year", "NAM", "ANG", "Total"))
  expect_equal(tab$NAM[tab$year == 1985], 2918)
  # close to the JARA (JAGS) reconstruction: ANG 1988 = 1923.4, NAM 2010 = 833.2
  expect_lt(abs(tab$ANG[tab$year == 1988]/1923.4 - 1), 0.1); expect_lt(abs(tab$NAM[tab$year == 2010]/833.2 - 1), 0.1)
  expect_equal(cf$catch, tab$Total)
  tr <- jb3_trend(jara_examples$monk[jara_examples$monk$year >= 1988, c("year", "survey1", "survey2")], timeblock = 2005)
  expect_equal(nrow(tr$r), 2)
})

test_that("Conn (2010) independent year effects", {
  data(iccat)
  a <- jb3_jointindex(jb3_trend(iccat$bum$cpue, iccat$bum$se, fixed.obsE = 0.1))
  tb <- jb3_trend(iccat$bum$cpue, iccat$bum$se, fixed.obsE = 0.1, state = "iid")
  expect_equal(tb$settings$state, "iid"); expect_equal(tb$convergence$tier, "ok")
  b <- jb3_jointindex(tb); ok <- is.finite(a$index) & is.finite(b$index)
  expect_gt(cor(log(a$index[ok]), log(b$index[ok])), 0.98)
})

test_that("catchfill with per-value SE and the DPS example", {
  data(jara_examples); x <- jara_examples$dps
  expect_true(all(c("dps", "monk", "dentex", "ank") %in% names(jara_examples)))
  cc <- x$catch; cc[1, "fleet5"] <- cc[1, "fleet2"]; cc[1, "fleet4"] <- 10
  sc <- cc; sc[, -1] <- 0.05; sc[1, c("fleet4", "fleet5")] <- 0.2
  cf <- jb3_catchfill(cc, se = sc, timeblock = 2016)
  expect_s3_class(cf, "jb3catchfill"); expect_equal(nrow(cf), nrow(cc))
  expect_gt(max(cf$se.log[cf$filled]), max(cf$se.log[!cf$filled]))
  cf0 <- jb3_catchfill(cc, se = 0.05, timeblock = 2016)               # scalar still works
  expect_true(all(is.finite(cf0$catch)))
  bi <- jb3_trend(x$bio, x$bio.se, type = "census")
  expect_equal(nrow(bi$total), nrow(x$bio))
})
