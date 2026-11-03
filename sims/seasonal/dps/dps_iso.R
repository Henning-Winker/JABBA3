`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs)
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
dS <- jb3_data(catch.s, data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs), data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE), nseason = 4)
prK <- jb3_priors(dS, r = c(1.07, 0.3), bkinit = c(0.9, 0.3)); prK$sdb <- c(median = 0.15, log.sd = 0.5, on = 1); prK$K[3] <- 0
for(ns in c(1, 4)) for(x in list(list(lag = 0, rho = 0), list(lag = 1, rho = 0), list(lag = 0, rho = 0.35), list(lag = 1, rho = 0.35),
                                 list(lag = 1, rho = 0.35, start = "index"), list(lag = 1, rho = 0, start = "index"))){
  lg <- if(ns == 4 && x$lag > 0) 0.9 else x$lag
  f <- tryCatch(suppressWarnings(jabba3(dS, prK, jb3_model(shape = "schaefer", nseason = ns, lag = lg, rho = x$rho, sdc = 0.2, start = x$start %||% "all"), osa = FALSE)), error = function(e) NULL)
  if(is.null(f)){ cat("error\n"); next }
  cat(sprintf("ns=%d lag %.1f rho %.2f start %-5s: %-7s maxgrad %.1e | K %6.0f r %.2f | B/Bmsy 2021 %.2f F/Fmsy %.2f\n", ns, lg, x$rho, x$start %||% "all",
      f$convergence$tier, f$convergence$maxgrad, f$report$K*f$mucatch, f$report$r, tail(f$report$BBmsy, 1), tail(f$report$HHmsy, 1)))
}
