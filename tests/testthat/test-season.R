test_that("nseason = 1 is the annual model; seasonal steps run and keep catch totals", {
  skip_if_not_installed("RTMB")
  data(iccat)
  d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
  pr <- jb3_priors(d, r = c(0.2, 0.5))
  f1 <- suppressWarnings(jabba3(d, pr, jb3_model(shape = "fox", lag = 2, rho = 0.5), osa = FALSE))
  f1b <- suppressWarnings(jabba3(d, pr, jb3_model(shape = "fox", lag = 2, rho = 0.5, nseason = 1), osa = FALSE))
  expect_equal(f1$opt$objective, f1b$opt$objective)
  expect_null(f1$ts.season)

  for(alloc in list("F", "equal", c(0.1, 0.4, 0.4, 0.1))){
    f4 <- suppressWarnings(jabba3(d, pr, jb3_model(shape = "fox", lag = 2, rho = 0.5, nseason = 4, catch.alloc = alloc), osa = FALSE))
    expect_equal(f4$nseason, 4L)
    expect_true(f4$convergence$tier %in% c("ok", "warning"))
    cs <- f4$ts.season[f4$ts.season$quant == "Catch", ]
    ca <- f4$ts[f4$ts$quant == "Catch", ]
    expect_equal(as.numeric(tapply(cs$est, cs$year, sum)), ca$est, tolerance = 1e-6)    # seasons add up to the annual catch
    expect_equal(nrow(f4$ts[f4$ts$quant == "B/Bmsy", ]), length(d$years))           # annual outputs
    expect_equal(f4$report$vf, .season_vfac(4, 0.5^(1/4)))
  }
  # legacy state-based equal F ("Fstate"): within-year harvest fraction constant on the states
  f4s <- suppressWarnings(jabba3(d, pr, jb3_model(shape = "fox", nseason = 4, catch.alloc = "Fstate"), osa = FALSE))
  h <- matrix(f4s$report$HS, 4)
  expect_lt(max(apply(h, 2, function(x) diff(range(x))/mean(x))), 1e-6)
  # sequential equal F (default): shares depend only on the start-of-year state; seasons add up to the annual catch
  f4 <- suppressWarnings(jabba3(d, pr, jb3_model(shape = "fox", nseason = 4), osa = FALSE))
  cs <- matrix(f4$report$estCS, 4)
  expect_equal(colSums(cs), as.numeric(f4$report$estC), tolerance = 1e-8)
  # annual model status close to the seasonal one for a slow stock
  bb1 <- tail(f1$ts$est[f1$ts$quant == "B/Bmsy"], 1); bb4 <- tail(f4$ts$est[f4$ts$quant == "B/Bmsy"], 1)
  expect_lt(abs(log(bb4/bb1)), 0.3)
})

test_that("seasonal inputs: season and time columns", {
  yrs <- 2001:2010
  cs <- expand.grid(season = 1:4, year = yrs)[, 2:1]; cs$catch <- 10
  cs <- rbind(data.frame(year = 2000, season = NA, catch = 40), cs)
  idx <- data.frame(time = rep(yrs, each = 2) + c(0.25, 0.75), survey = exp(rnorm(20, 0, 0.1)))
  d <- jb3_data(cs, idx, nseason = 4, stock = "sea")
  expect_equal(d$nseason, 4L)
  expect_equal(d$catch, rep(40, 11))
  expect_true(all(is.na(d$catch.season[1, ])))
  expect_equal(nrow(d$iobs), 20)
  st <- .season_steps(d, jb3_model(nseason = 4), d$years)
  expect_equal(sum(!is.na(st$I)), 20)
  expect_equal(unique(st$stepseason[!is.na(st$I[, 1])]), c(2, 4))                     # 0.25 -> season 2, 0.75 -> season 4
  expect_equal(st$allocF, c(1L, rep(0L, 10)))                                     # seasonal catch given from 2001
  expect_error(jb3_data(cs[-3, ], idx, nseason = 4), "incomplete")
  # annual index with timing
  d2 <- jb3_data(data.frame(year = yrs, catch = 40), data.frame(year = yrs, cpue = 1), index.time = 0.5)
  st2 <- .season_steps(d2, jb3_model(nseason = 4), d2$years)
  expect_equal(unique(st2$stepseason[!is.na(st2$I[, 1])]), 3)
  expect_error(jb3_model(nseason = 4, season.w = c(1, 2)), "season.w")
})
