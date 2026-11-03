`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs)
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
dS <- jb3_data(catch.s, data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs), data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE), nseason = 4)
pr <- jb3_priors(dS, r = c(1.072792*exp(-0.5*0.3^2), 0.3), bkinit = c(0.9*exp(-0.5*0.3^2), 0.3)); pr$sdb <- c(median = 0.15*exp(-0.5*0.5^2), log.sd = 0.5, on = 1)
pr$K[3] <- 0                                                     # K free
for(ns in c(1, 4)) for(cap in c(NA, 0.3, 0.22)){
  md <- jb3_model(shape = "schaefer", nseason = ns, sdc = 0.2, sdi.est = FALSE, sdb.max = if(is.na(cap)) NULL else cap)
  f <- tryCatch(suppressWarnings(jabba3(dS, pr, md, osa = FALSE, nstart = 6)), error = function(e) NULL)
  if(is.null(f)){ cat("error\n"); next }
  bb <- f$ts[f$ts$quantity == "B/Bmsy", ]; ff <- f$ts[f$ts$quantity == "F/Fmsy", ]; kk <- f$refpts["K", ]; n <- nrow(bb)
  cat(sprintf("%-9s cap %-4s: %-4s K %6.0f (%6.0f-%6.0f) r %.2f sdb %.3f | B/Bmsy 2021 %.2f (%.2f-%.2f) F/Fmsy %.2f (%.2f-%.2f)\n", if(ns == 1) "annual" else "quarterly", cap,
      f$convergence$tier, kk$est, kk$lci, kk$uci, f$report$r, f$report$sdb, bb$est[n], bb$lci[n], bb$uci[n], ff$est[n], ff$lci[n], ff$uci[n]))
}
# likelihood profile over K (quarterly, fixed SE; no cap and cap 0.22)
Ks <- c(20e3, 25e3, 30e3, 35e3, 40e3, 45e3, 50e3, 60e3, 80e3, 120e3, 200e3)
for(cap in c(NA, 0.22)){
  pp <- sapply(Ks, function(K){ p <- pr; p$K <- c(K, 0.001, 1)
    md <- jb3_model(shape = "schaefer", nseason = 4, sdc = 0.2, sdi.est = FALSE, sdb.max = if(is.na(cap)) NULL else cap)
    f <- tryCatch(suppressWarnings(jabba3(dS, p, md, osa = FALSE, se = FALSE)), error = function(e) NULL)
    if(is.null(f)) return(c(NA, NA, NA)); c(f$opt$objective + dnorm(log(f$report$K*f$mucatch), log(K), 0.001, log = TRUE), tail(f$report$BBmsy, 1), f$report$r) })
  cat("\nprofile quarterly fixed SE, cap", cap, "\n"); print(round(rbind(K_thousands = Ks/1000, dNLL = pp[1, ] - min(pp[1, ], na.rm = TRUE), BBmsy2021 = pp[2, ], r = pp[3, ]), 2))
}
