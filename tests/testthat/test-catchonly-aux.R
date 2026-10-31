data(iccat, package = "JABBA3")
bet <- iccat$bet

test_that("catch-only data, status priors and the default NUTS engine", {
  d0 <- jb3_data(bet$catch, NULL, stock = "bet")
  expect_equal(ncol(d0$I), 0)
  expect_output(print(d0), "catch-only")
  st <- jb3_bprior(0.7, 0.2)
  expect_true(is.na(st$year)); expect_equal(st$type, "BBmsy"); expect_equal(st$log.sd, sqrt(log(1 + 0.2^2)))
  pr <- jb3_priors(d0, r = c(0.2, 0.5), status = c(0.5, 0.2, 2010, "bbmsy"))     # JABBA b.prior vector
  expect_equal(pr$status$year, 2010); expect_equal(pr$status$log.sd, 0.2)
  pr0 <- jb3_priors(d0, r = c(0.2, 0.5), status = rbind(st, jb3_bprior(1.2, 0.3, 2000, "ffmsy"), jb3_bprior(0.6, 0.2, 1990, "bk")))
  s <- JABBA3:::jb3_setup(d0, pr0, jb3_model(shape = "fox"))
  expect_equal(s$dat$st.t, c(68, 51, 41))
  expect_equal(s$dat$nI, 0); expect_true(all(is.na(s$map$logq)))
  f <- suppressWarnings(jabba3(d0, pr0, jb3_model(shape = "fox"), engine = "rtmb", se = FALSE))
  expect_equal(nrow(f$index), 0); expect_false("q1" %in% rownames(f$pars))
  expect_message(fn <- suppressWarnings(jabba3(d0, pr0, jb3_model(shape = "fox"), nuts = list(chains = 1, iter = 200))), "nuts")
  skip_if_not_installed("tmbstan")
  expect_equal(fn$engine, "nuts"); expect_equal(fn$nuts$param, "centred")
})

test_that("non-centred states reproduce the centred process model", {
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  for(m in list(jb3_model(shape = "fox"), jb3_model(shape = "fox", rho = 0.6, lag = 2, biascor = TRUE),
                jb3_model(shape = "fox", rw = TRUE, start = "index"), jb3_model(shape = "fox", sdb.form = "dyn"))){
    f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), m, osa = FALSE, se = FALSE)
    o <- JABBA3:::.nc_obj(f)
    expect_equal(o$report(o$env$par)$P, f$report$P, tolerance = 1e-8)
  }
  # mode fits are unchanged by the parameterisation (always centred)
  f1 <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox", param = "noncentred"), osa = FALSE)
  expect_equal(f1$pars["K", "est"], 1174749, tolerance = 1e-4)
})

test_that("prior-only runs switch off the data likelihoods", {
  d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
  s <- JABBA3:::jb3_setup(d, jb3_priors(d), jb3_model(), prior.only = TRUE)
  expect_true(s$dat$prior.only); expect_true(all(is.na(s$map$logq)))
  obj <- RTMB::MakeADFun(function(p) JABBA3:::jb3_nll(p, s$dat), s$par, map = s$map, silent = TRUE)
  d2 <- d; d2$I[] <- d2$I*runif(length(d2$I), 0.5, 2)                       # different index values, same nll
  s2 <- JABBA3:::jb3_setup(d2, jb3_priors(d), jb3_model(), prior.only = TRUE)
  obj2 <- RTMB::MakeADFun(function(p) JABBA3:::jb3_nll(p, s2$dat), s$par, map = s$map, silent = TRUE)
  expect_equal(obj$fn(obj$par), obj2$fn(obj$par))
})

test_that("auxiliary data: input formats, effort creep and fits", {
  expect_true(all(c("effort_hyp", "ffmsy_hyp") %in% names(bet)))
  da <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet", auxiliary = bet$effort_hyp, auxiliary.type = "effort")
  expect_equal(unname(da$A.type), "effort"); expect_equal(sum(!is.na(da$A)), nrow(bet$effort_hyp))
  dl <- jb3_data(bet$catch, NULL, auxiliary = list(effort = bet$effort_hyp, ffmsy = list(data = bet$ffmsy_hyp[, 1:2], se = 0.3, type = "ffmsy")))
  expect_equal(unname(dl$A.type), c("effort", "ffmsy")); expect_equal(unname(dl$A.SE[as.character(1995), "ffmsy"]), 0.3)
  expect_error(jb3_data(bet$catch, NULL, auxiliary = bet$effort_hyp, auxiliary.type = "cpue"), "auxiliary.type")
  pr <- jb3_priors(da, r = c(0.2, 0.5))
  f0 <- jabba3(da, pr, jb3_model(shape = "fox", effort.creep = 0.02), osa = FALSE)
  expect_equal(f0$convergence$tier, "ok")
  expect_true(all(c("qA.effort", "sdA.effort") %in% rownames(f0$pars)))
  expect_equal(nrow(f0$auxiliary), nrow(bet$effort_hyp))
  fe <- jabba3(da, pr, jb3_model(shape = "fox", effort.creep = "estimate"), osa = FALSE)
  expect_true("creep" %in% rownames(fe$pars))
  expect_lt(abs(fe$pars["creep", "est"] - 0.02), 0.02)       # simulated with 2% creep
  ff <- jabba3(da, pr, jb3_model(shape = "fox", effort.creep = 0.02, aux.sd = 0.15), osa = FALSE)
  expect_equal(ff$pars["sdA.effort", "est"], 0.15)
  expect_s3_class(jbplot_auxfits(f0), "ggplot")
  # status-type auxiliary data (F/Fmsy), q fixed at 1
  dF <- jb3_data(bet$catch, NULL, auxiliary = bet$ffmsy_hyp[, 1:2], auxiliary.se = bet$ffmsy_hyp[, c(1, 3)], auxiliary.type = "ffmsy")
  s <- JABBA3:::jb3_setup(dF, jb3_priors(dF), jb3_model(shape = "fox"))
  expect_true(is.na(s$map$logqA)); expect_equal(s$dat$atype, 3L)
})

test_that("JABBA converters keep auxiliary data and catch-only", {
  skip_if_not_installed("FLCore")
  x <- jbdat2flr(bet$catch, bet$cpue, bet$se, auxiliary = bet$effort_hyp)
  expect_s4_class(x$auxiliary, "FLQuants")
  b <- flr2jbdat(x$catch, x$indices, x$auxiliary)
  expect_equal(b$auxiliary$effort[b$auxiliary$year == 1990], bet$effort_hyp$effort[bet$effort_hyp$year == 1990])
  d <- jb3_data(x$catch, x$indices, stock = "bet", auxiliary = x$auxiliary)
  expect_equal(colnames(d$A), "effort")
  c0 <- flr2jbdat(x$catch, NULL, x$auxiliary); expect_null(c0$cpue)
})

test_that("jbinput2jabba3 maps auxiliary and catch-only JABBA inputs", {
  skip_if_not_installed("JABBA")
  jb <- suppressMessages(JABBA::build_jabba(catch = bet$catch, cpue = bet$cpue, se = bet$se, model.type = "Fox",
    auxiliary = merge(bet$catch[, 1, drop = FALSE], bet$effort_hyp, all.x = TRUE), auxiliary.type = "effort", verbose = FALSE))
  x <- suppressMessages(jbinput2jabba3(jb, verbose = FALSE))
  expect_equal(unname(x$data$A.type), "effort")
  jc <- suppressMessages(JABBA::build_jabba(catch = bet$catch, model.type = "Fox", scenario = "CatchOnly",
    b.prior = c(0.5, 0.2, 2010, "bbmsy"), psi.prior = c(1, 0.1), verbose = FALSE))
  y <- suppressMessages(jbinput2jabba3(jc, verbose = FALSE))
  expect_equal(ncol(y$data$I), 0)
  expect_equal(y$priors$status$type, "BBmsy"); expect_equal(y$priors$status$year, 2010)
})
