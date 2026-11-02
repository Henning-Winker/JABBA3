`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.a <- aggregate(Obs ~ Yr, cs, sum); names(catch.a) <- c("year", "catch")
prof <- tapply(cs$Obs, cs$Seas, sum); prof <- as.numeric(prof/sum(prof)); cat("SS3 mean seasonal shares:", round(prof, 3), "\n")
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean)
dA <- jb3_data(catch.a, data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs), data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE), nseason = 4)
pr <- jb3_priors(dA, r = c(1.07, 0.3), bkinit = c(0.9, 0.3)); pr$sdb <- c(median = 0.15, log.sd = 0.5, on = 1)
for(al in list("F", "Fstate", "equal", prof)){
  f <- suppressWarnings(jabba3(dA, pr, jb3_model(shape = "schaefer", nseason = 4, sdc = 0.2, catch.alloc = al), osa = FALSE))
  cat(sprintf("%-28s K %6.0f  r %.2f  B/Bmsy 2021 %.2f  F/Fmsy 2021 %.2f  nll %.2f  %s\n", paste(if(is.numeric(al)) "shares SS3 profile" else al), f$report$K*f$mucatch, f$report$r,
      tail(f$report$BBmsy, 1), tail(f$report$HHmsy, 1), f$opt$objective, f$convergence$tier)) }
