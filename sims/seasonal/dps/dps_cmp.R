`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs)
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
dS <- jb3_data(catch.s, data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs), data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE), nseason = 4)
oe <- mean(ci$SE)
mk <- function(spict){ p <- jb3_priors(dS, r = c(1.072792*exp(-0.5*0.3^2), 0.3), bkinit = c(0.9*exp(-0.5*0.3^2), 0.3)); p$sdb <- c(median = 0.15*exp(-0.5*0.5^2), log.sd = 0.5, on = 1); p$K[3] <- 0
  if(spict) p$sdi <- c(median = oe, log.sd = 0.3, on = 1); p }
sp <- readRDS("spict_ref.rds")
out <- list()
for(ns in c(1, 4)) for(obs in c("additive", "spict")){
  f <- suppressWarnings(jabba3(dS, mk(obs == "spict"), jb3_model(shape = "schaefer", nseason = ns, sdc = 0.2, obs.error = obs), osa = FALSE))
  bb <- f$ts[f$ts$quantity == "B/Bmsy", ]; ff <- f$ts[f$ts$quantity == "F/Fmsy", ]; n <- nrow(bb)
  cat(sprintf("JABBA3 ns=%d obs %-8s: %-4s K %6.0f r %.2f sdb %.3f sdi %.3f | 2021 B/Bmsy %.2f (%.2f-%.2f) F/Fmsy %.2f (%.2f-%.2f)\n", ns, obs, f$convergence$tier,
      f$report$K*f$mucatch, f$report$r, f$report$sdb, f$report$sdi[1], bb$est[n], bb$lci[n], bb$uci[n], ff$est[n], ff$lci[n], ff$uci[n]))
  out[[paste(ns, obs)]] <- data.frame(year = bb$year, BB = bb$est, FF = ff$est)
}
cat("SPiCT                    :      K", round(sp$K), " r", round(sp$r, 2), "\n")
yrs <- c(1970, 1980, 1990, 2000, 2005, 2010, 2015, 2018, 2020, 2021)
spB <- sapply(yrs, function(y) sp$bb[which.min(abs(sp$time - y)), 2]); spF <- sapply(yrs, function(y) sp$ff[which.min(abs(sp$time - (y + 0.5))), 2])
tb <- data.frame(year = yrs, SPiCT = spB, sapply(out, function(o) o$BB[match(yrs, o$year)]), check.names = FALSE)
tf <- data.frame(year = yrs, SPiCT = spF, sapply(out, function(o) o$FF[match(yrs, o$year)]), check.names = FALSE)
cat("\nB/Bmsy (JABBA3 start of year; SPiCT at Jan 1)\n"); print(tb, digits = 2, row.names = FALSE)
cat("\nF/Fmsy (SPiCT mid-year instantaneous; JABBA3 annual harvest rate / Hmsy)\n"); print(tf, digits = 2, row.names = FALSE)
