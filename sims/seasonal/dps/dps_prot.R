`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs)
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
dS <- jb3_data(catch.s, data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs), data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE), nseason = 4)
prK <- jb3_priors(dS, r = c(1.07, 0.3), bkinit = c(0.9, 0.3)); prK$sdb <- c(median = 0.15, log.sd = 0.5, on = 1); prK$K[3] <- 0
for(x in list()){
  f <- suppressWarnings(jabba3(dS, prK, jb3_model(shape = "schaefer", nseason = x$ns, lag = x$lag, rho = 0.35, sdc = 0.2, sdb.max = x$cap), osa = FALSE, nstart = 8))
  cv <- f$convergence
  cat(sprintf("ns=%d lag %.1f cap %-4s: tier %-7s pdHess %s maxgrad %.1e | K %6.0f r %.2f sdb %.3f | B/Bmsy 2021 %.2f F/Fmsy %.2f\n", x$ns, x$lag, x$cap %||% "-", cv$tier,
      cv$pdHess, cv$maxgrad, f$report$K*f$mucatch, f$report$r, f$report$sdb, tail(f$report$BBmsy, 1), tail(f$report$HHmsy, 1)))
  if(!isTRUE(cv$pdHess)){ h <- tryCatch(optimHess(f$opt$par, f$obj$fn, f$obj$gr), error = function(e) NULL)
    if(!is.null(h)){ ev <- eigen(h, symmetric = TRUE); k <- which.min(ev$values); v <- ev$vectors[, k]; names(v) <- names(f$opt$par)
      cat("   smallest Hessian eigenvalue", signif(ev$values[k], 3), "; direction:", paste(names(v)[order(-abs(v))][1:3], round(v[order(-abs(v))][1:3], 2), collapse = ", "), "\n") } }
}
cat("\n--- where is it non-smooth? ---\n")
for(x in list(list(ns = 1, lag = 1), list(ns = 4, lag = 0.9))){
  f <- suppressWarnings(jabba3(dS, prK, jb3_model(shape = "schaefer", nseason = x$ns, lag = x$lag, rho = 0.35, sdc = 0.2), osa = FALSE))
  rp <- f$report; d <- f$setup$dat
  P <- rp$PS; Pl <- P[d$lagidx[-d$N]]
  arg <- P[-d$N] + d$gfac[-d$N]*rp$r/(rp$m - 1)*Pl*(1 - Pl^(rp$m - 1)) - rp$estCS[-d$N]/(rp$K*rp$remS[-d$N])
  cat(sprintf("ns=%d: P range %.3f-%.3f (bounds %.2f-%.2f); steps with P > 1.3: %d; production argument min %.4f (floor 0.001), steps <= 0.002: %d\n",
      x$ns, min(P), max(P), d$P.bound[1], d$P.bound[2], sum(P > d$P.bound[2]), min(arg), sum(arg <= 0.002)))
  w <- which(P > 1.25); if(length(w)) cat("   steps near/above upper bound:", head(w, 12), " (first years:", head(unique(f$years[d$stepyear[w]]), 6), ")\n")
}
