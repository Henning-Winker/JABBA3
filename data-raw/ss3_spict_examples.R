# Builds data/ss3_examples.rda and data/spict_examples.rda from the BioProcDyn analyses (branch wkbmsyspict):
#   ss3_examples: North Sea sole (sol.27.4) and blackspot seabream (sbr, GFCM GSA 1 and 3): JABBA3 inputs (data frames),
#                 priors used in the guide, and the SS3 reference trajectories (data.frame(year, quant, est, lci, uci))
#   spict_examples: WKBMSYSPiCT3 benchmark inputs (SPiCT inp lists) for nep FU32, aru.27.1234a, dab.27.22-32 and
#                 cod.27.2.coastS, with the SPiCT fit trajectories in the same format
# Run from the JABBA3 root with BioProcDyn at ../bioprocdyn: Rscript data-raw/ss3_spict_examples.R
suppressMessages({library(FLCore); library(spict); library(data.table)})
bp <- "../bioprocdyn"
## North Sea sole
fl <- readRDS(file.path(bp, "ss3sol/data/sol.27.4_flr.rds"))
fq <- function(x) data.frame(year = as.numeric(dimnames(x)$year), value = c(x))
ct <- fq(fl$catch.total); names(ct)[2] <- "catch"
ix <- Reduce(function(a, b) merge(a, b, all = TRUE), lapply(names(fl$indices), function(n){ x <- fq(index(fl$indices[[n]])); names(x)[2] <- n; x }))
sx <- Reduce(function(a, b) merge(a, b, all = TRUE), lapply(names(fl$indices), function(n){ x <- fq(index.var(fl$indices[[n]])); names(x)[2] <- n; x }))
s3 <- as.data.frame(readRDS(file.path(bp, "ss3sol/output/sol_runs.rds"))$ss3)[, c("year", "quant", "est", "lci", "uci")]
sol <- list(stock = "North Sea sole (sol.27.4)", catch = ct, cpue = ix, se = sx,
  info = list(assessment = "ICES SS3 3.30.22 (2025): single sex, ages 0-10+, empirical weight-at-age, landings + discard fleets, BTS and COAST surveys (numbers) 1990-2025",
              catch = "total dead catch (landings + dead discards), t", indices = "survey biomass indices from the SS3 data (ss3diags2::ss2jabba3), log-scale SEs",
              priors = list(r = c(0.181, 0.206), BmsyK = 0.336, lag = 2, rho = 0.72, note = "age-structured equilibrium of the 1970-80 benchmark biology (EB state); GT 6.6, a_mat 2"),
              refpts = c(SSB0 = 276000, SSBmsy = 67100, Fmsy = 0.301, MSY = 20000), h = fl$h),
  biology = fl$bio.at.age, ss3 = s3)
## blackspot seabream
sb <- readRDS(file.path(bp, "ss3sbr/example/sbr_jabba3_input.rds"))
s3b <- as.data.frame(sb$ss3)[, c("year", "quant", "est", "lci", "uci")]
sbr <- list(stock = "Blackspot seabream (Pagellus bogaraveo), GFCM GSA 1 and 3 (sbr.37.1-3)", catch = sb$catch, cpue = sb$cpue, se = sb$se,
  info = list(assessment = "GFCM SS3 3.30.22 (2025 update, 1983-2024): two sexes, protandric hermaphrodite, M-at-age, h 0.7, ESP and MOR fleets and CPUE (2009-2024)",
              catch = "total dead catch of both fleets, t", indices = "ESP and MOR fishery CPUE with input SEs; MOR down-weighted (SE x 1.5) in the guide (potential hyperstability)",
              priors = list(r = c(0.157, 0.22), BmsyK = 0.4, lag = 4, rho = 0.82, note = "r = m x H_MSY from the age-structured equilibrium of the SS3 biology; lag = age of 50% mature females; GT 11.3"),
              refpts = sb$ss3_refpts, forecast = sb$ss3_forecast, h = sb$h),
  biology = sb$biology, ss3 = s3b)
ss3_examples <- list(sol = sol, sbr = sbr)
## SPiCT stocks
files <- c(nep = "nep.fu.32_inp_reconstructed.rds", aru = "aru.27.1234a_inp_final.rds", dab = "dab.27.22-32_inp_final.rds", cod = "cod.27.2.coastS_inp_final.rds")
names_full <- c(nep = "nep FU32", aru = "aru.27.1234a", dab = "dab.27.22-32", cod = "cod.27.2.coastS")
ann <- function(fit, par){ x <- get.par(par, fit, exp = TRUE); t <- as.numeric(rownames(x))
  d <- data.table(time = t, est = x[, "est"], lci = x[, "ll"], uci = x[, "ul"]); d[, year := floor(time + 1e-8)]
  if(par == "logBBmsy") d[abs(time - year) < 1e-6, .(year, est, lci, uci)] else d[, .(est = mean(est), lci = mean(lci), uci = mean(uci)), by = year] }
spict_examples <- lapply(names(files), function(s){ inp <- readRDS(file.path(bp, "wkbmsyspict/data/wk3", files[[s]])); meta <- inp$meta; inp$meta <- NULL
  fit <- fit.spict(inp); yl <- floor(max(inp$timeC))
  ts <- rbind(data.table(quant = "B/Bmsy", ann(fit, "logBBmsy")), data.table(quant = "F/Fmsy", ann(fit, "logFFmsy")))[year <= yl, .(year, quant, est, lci, uci)]
  list(stock = names_full[[s]], inp = inp, spict = as.data.frame(ts),
       info = list(benchmark = "ICES WKBMSYSPiCT3 (2024), accepted run", K = get.par("logK", fit, exp = TRUE)[, "est"], r = get.par("logr", fit, exp = TRUE)[, "est"],
                   maxC = max(inp$obsC), meta = meta)) })
names(spict_examples) <- names(files)
save(ss3_examples, file = "data/ss3_examples.rda", compress = "xz")
save(spict_examples, file = "data/spict_examples.rda", compress = "xz")
