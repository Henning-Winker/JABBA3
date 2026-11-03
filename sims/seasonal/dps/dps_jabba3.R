suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
`%||%` <- function(a, b) if(is.null(a)) b else a
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
# seasonal catch (sum over fleets), season from SS3
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs)
catch.a <- aggregate(catch ~ year, catch.s, sum)
# survey: SS3 season, mean of duplicates as in the SPiCT Rmd
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
idx.s <- data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs); se.s <- data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE)
cat("catch years", range(catch.a$year), " survey obs", nrow(idx.s), "\n")
dS <- jb3_data(catch.s, idx.s, se.s, nseason = 4, stock = "DPS 12-16")           # seasonal catch + seasonal survey
dA <- jb3_data(catch.a, idx.s, se.s, nseason = 4, stock = "DPS 12-16")           # annual catch + seasonal survey
pr <- function(d){ p <- jb3_priors(d, r = c(1.07, 0.3), bkinit = c(0.9, 0.3)); p$sdb <- c(median = 0.15, log.sd = 0.5, on = 1); p }
# protocol for a short-lived stock: lag at a50 maturity (0.9 y), rho 0.35 (GT <= 3)
runs <- list(
  annual            = list(d = dA, ns = 1),
  annual_prot       = list(d = dA, ns = 1, lag = 1, rho = 0.35),
  q_seascatch       = list(d = dS, ns = 4),
  q_seascatch_prot  = list(d = dS, ns = 4, lag = 0.9, rho = 0.35),
  q_seascatch_prot_cap = list(d = dS, ns = 4, lag = 0.9, rho = 0.35, cap = 0.3),
  q_anncatch        = list(d = dA, ns = 4),
  q_anncatch_prot   = list(d = dA, ns = 4, lag = 0.9, rho = 0.35))
fits <- lapply(runs, function(x){
  md <- jb3_model(shape = "schaefer", nseason = x$ns, lag = x$lag %||% 0, rho = x$rho %||% 0, sdc = 0.2, sdb.max = x$cap)
  tryCatch(suppressWarnings(jabba3(x$d, pr(x$d), md, osa = FALSE)), error = function(e){ message(conditionMessage(e)); NULL }) })
saveRDS(fits, "dps_fits.rds")
# SS3 status
d <- r$derived_quants; yrs <- 1970:2021
ssb <- d[paste0("SSB_", yrs), "Value"]; ssbmsy <- d["SSB_MSY", "Value"]; Fv <- d[paste0("F_", yrs), "Value"]; fmsy <- d["annF_MSY", "Value"]
ss3 <- data.frame(year = yrs, BBmsy = ssb/ssbmsy, FFmsy = Fv/fmsy)
saveRDS(ss3, "ss3_status.rds")
tab <- do.call(rbind, lapply(names(fits), function(k){ f <- fits[[k]]; if(is.null(f)) return(data.frame(run = k, conv = "error"))
  g <- function(q) f$ts[f$ts$quantity == q, ]; bb <- g("B/Bmsy"); ff <- g("F/Fmsy"); n <- nrow(bb)
  data.frame(run = k, conv = f$convergence$tier, r = round(f$report$r, 2), sdb = round(f$report$sdb, 3), Fmsy = round(f$report$Hmsy, 2),
             BBmsy2021 = sprintf("%.2f (%.2f-%.2f)", bb$est[n], bb$lci[n], bb$uci[n]), FFmsy2021 = sprintf("%.2f (%.2f-%.2f)", ff$est[n], ff$lci[n], ff$uci[n]),
             BBmsy_min = round(min(bb$est), 2), nll = round(f$opt$objective, 1)) }))
print(tab, row.names = FALSE)
cat("\nSS3 (SSB/SSBmsy, F/Fmsy) 2021:", round(tail(ss3$BBmsy, 1), 2), round(tail(ss3$FFmsy, 1), 2), " min SSB/SSBmsy", round(min(ss3$BBmsy), 2), "\n")
cat("SPiCT reference (from the trial report): B/Bmsy ~0.96 (2023 forecast, keep current catch), F/Fmsy ~1.0\n")
