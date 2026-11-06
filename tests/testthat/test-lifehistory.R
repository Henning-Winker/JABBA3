test_that("Leslie r and GT reproduce JABBA::jbleslie", {
  skip_if_not_installed("JABBA")
  bio <- list(age = 0:15, M = 0.2, Mat = 1/(1 + exp(-(0:15 - 3))), Wa = (1 - exp(-0.2*((0:15) + 0.5)))^3, h = 0.7)
  j <- JABBA::jbleslie(bio)
  m <- jb3_leslie(list(age = bio$age, M = bio$M, mat = bio$Mat, wt = bio$Wa, h = 0.7))
  expect_equal(unname(m["r"]), j$r, tolerance = 1e-10); expect_equal(unname(m["gt"]), j$gt, tolerance = 1e-10)
})

test_that("Leslie on FLStock reproduces FLSRTMB::productivity (constant M)", {
  skip_if_not_installed("FLSRTMB"); skip_if_not_installed("FLCore")
  data(ple4, package = "FLCore", envir = environment())
  s <- ple4[1:10, "2000"]; FLCore::m(s)[] <- 0.2; FLCore::m.spwn(s)[] <- 0
  pf <- FLSRTMB::productivity(s, s = 0.7)
  m <- jb3_leslie(s, h = 0.7)
  expect_equal(unname(m["r"]), c(pf$r), tolerance = 1e-8); expect_equal(unname(m["gt"]), c(pf$gt), tolerance = 1e-8)
})

test_that("lhpar, sex structure and Fmsy-r conversions", {
  skip_if_not_installed("FLCore")
  lh <- FLCore::FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.2, a = 0.0047, b = 3.15, s = 0.75)
  r1 <- jb3_leslie(lh); expect_true(r1["r"] > 0 && r1["gt"] > 1)
  expect_equal(jb3_leslie(c(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.2, a = 0.0047, b = 3.15, s = 0.75)), r1)
  lh2 <- FLCore::FLPar(array(c(55.7, 0.08, 0.162, 31.2, 0.0047, 3.15, 0.75, 50, 0.09, 0.18, 29, 0.0047, 3.15, 0.75),
    dim = c(7, 2, 1), dimnames = list(params = c("linf", "k", "M", "L50", "a", "b", "s"), sex = c("F", "M"), iter = 1)))
  expect_equal(jb3_leslie(lh2, sex = "F")[["r"]], jb3_leslie(lh)[["r"]], tolerance = 1e-3)
  expect_false(isTRUE(all.equal(jb3_leslie(lh2, sex = "M"), jb3_leslie(lh2, sex = "F"))))
  expect_equal(jb3_fmsy2r(0.2, shape = "schaefer"), 0.4)
  expect_equal(jb3_r2fmsy(jb3_fmsy2r(0.2, BmsyK = 0.36), BmsyK = 0.36), 0.2)
  skip_if_not_installed("JABBA")
  expect_equal(jb3_fmsy2r(0.2, BmsyK = 0.36), JABBA::r_from_fmsy_bmsyk(0.2, 0.36)$r, tolerance = 1e-6)
})

test_that("Selex curves from FLicc length selectivity and the JABBA 5-parameter form", {
  lh <- c(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.2, a = 0.0047, b = 3.15, s = 0.75)
  sx <- jb3_selex(lh, list(G = c(mode = 32, lsd = 4, rsd = 15), L = c(SL50 = 22, SL95 = 28),
                           J = c(SL50 = 25, SL95 = 32, SLdesc = 45, CVdesc = 0.15, Smin = 0.1)))
  expect_s3_class(sx, "jb3selex"); expect_equal(names(sx$ebsb), c("G", "L", "J"))
  expect_true(all(sapply(sx$ebsb, length) == 5))
  cr <- sx$curves[sx$curves$P >= 0.1 & sx$curves$P <= 0.9, ]; expect_lt(max(abs(cr$fit - cr$EB_SB)), 0.15)   # fitted over 0.1-0.9
  sa <- jb3_selex(lh, list(A = c(SL50 = 3, SL95 = 5, SLdesc = 15, CVdesc = 0.2, Smin = 0.2)), sel.type = "age")
  expect_equal(length(sa$ebsb$A), 5)
  expect_error(jb3_selex(lh, list(c(SL50 = 3, SL95 = 5)), sel.type = "age"), "5-parameter")
  s5 <- jb3_sel5(c(0, 35, 45, 70, 200), c(SL50 = 35, SL95 = 45, SLdesc = 70, CVdesc = 0.15, Smin = 0.05))
  expect_equal(s5[2], 0.5); expect_equal(s5[4], 1, tolerance = 1e-3); expect_lt(s5[5], 0.1); expect_true(all(s5 >= 0 & s5 <= 1))
  expect_s3_class(jbplot_selex(sx), "patchwork"); expect_s3_class(jbplot_selex(sa, at = "age"), "patchwork")
  m <- jb3_model(ebsb = sx); expect_equal(length(m$ebsb), 3)
})

test_that("jb3_selex numbers curves and jb3_asem_prior", {
  lh <- c(linf = 147.5, k = 0.126, t0 = -1.89, a = 1.3718e-5, b = 3.0793, M = 0.3, s = 0.64, L50 = 79.9, L95 = 91.3, amax = 15)
  mat <- c(SL50 = 79.9, SL95 = 91.3); s1 <- c(SL50 = 57.2, SL95 = 72.4)
  sx <- jb3_selex(lh, list(LL = s1, LL.N = s1), units = c("biomass", "numbers"))
  expect_equal(unname(sx$units), c("biomass", "numbers"))
  P <- c(0.1, 0.9); w <- JABBA3:::.rho_curve(P, sx$ebsb$LL); n <- JABBA3:::.rho_curve(P, sx$ebsb$LL.N)
  expect_gt(n[1]/n[2], w[1]/w[2])                                      # numbers ratio rises more steeply with depletion
  d <- sx$curves[sx$curves$P >= 0.1 & sx$curves$P <= 0.9, ]; expect_lt(max(abs(d$fit/d$EB_SB - 1)), 0.05)
  expect_error(jb3_selex(lh, s1, units = "kg"), "units")
  expect_s3_class(jbplot_selex(sx), "patchwork")
  m <- jb3_model(ebsb = sx, ebsb.index = "LL.N", ebsb.catch = "LL"); expect_equal(length(m$ebsb), 2)
  # equilibrium prior: selection at maturity, catch removed from SB (H = MSY/SBmsy)
  p0 <- jb3_asem_prior(lh, mat, removal = "SB")
  expect_equal(nrow(p0$draws), 1); expect_equal(p0$BmsyK, 0.317, tolerance = 0.02); expect_equal(p0$draws$Hmsy, 0.238, tolerance = 0.02)
  expect_equal(unname(p0$r[1]), p0$draws$Hmsy*JABBA3:::.m_from_bmsyk(p0$BmsyK))
  p1 <- jb3_asem_prior(lh, s1, M.cv = 0.3, h.sd = 0.06, nsim = 50)
  expect_equal(nrow(p1$draws), 50); expect_gt(p1$r[2], 0); expect_s3_class(p1, "jb3asem")
  expect_equal(jb3_asem_prior(lh, s1, state = "EB")$draws$Hmsy, jb3_asem_prior(lh, s1)$draws$Hmsy)
  expect_output(print(p1), "r prior")
})
