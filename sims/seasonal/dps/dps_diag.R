`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs); catch.a <- aggregate(catch ~ year, catch.s, sum)
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
dS <- jb3_data(catch.s, data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs), data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE), nseason = 4)
pr0 <- jb3_priors(dS, r = c(1.07, 0.3), bkinit = c(0.9, 0.3)); pr0$sdb <- c(median = 0.15, log.sd = 0.5, on = 1)
cat("default K prior:", round(pr0$K, 2), "  max catch", max(catch.a$catch), "\n")
prK <- pr0; prK$K[3] <- 0
res <- list()
for(ns in c(1, 4)) for(kp in c("on", "off")){
  p <- if(kp == "on") pr0 else prK
  f <- suppressWarnings(jabba3(dS, p, jb3_model(shape = "schaefer", nseason = ns, sdc = 0.2), osa = FALSE))
  bb <- f$ts[f$ts$quantity == "B/Bmsy", ]; ff <- f$ts[f$ts$quantity == "F/Fmsy", ]; n <- nrow(bb)
  cat(sprintf("ns=%d K prior %-3s: K %6.0f r %.2f sdb %.3f sdi %.3f q %.2e | B/Bmsy 2021 %.2f (%.2f-%.2f) F/Fmsy %.2f (%.2f-%.2f) | %s\n", ns, kp,
      f$report$K*f$mucatch, f$report$r, f$report$sdb, f$report$sdi[1], f$report$q[1]/f$mucatch, bb$est[n], bb$lci[n], bb$uci[n], ff$est[n], ff$lci[n], ff$uci[n], f$convergence$tier))
  res[[paste(ns, kp)]] <- f
}
saveRDS(res, "dps_diag_fits.rds")
# likelihood profile over K (prior off), annual and quarterly: fix logK via map
prof <- function(ns, Ks){
  sapply(Ks, function(K){ p <- prK; p$K <- c(K, 0.001, 1)    # near-fixed K
    f <- tryCatch(suppressWarnings(jabba3(dS, p, jb3_model(shape = "schaefer", nseason = ns, sdc = 0.2), osa = FALSE, se = FALSE)), error = function(e) NULL)
    if(is.null(f)) return(c(NA, NA, NA)); c(f$opt$objective + dnorm(log(f$report$K*f$mucatch), log(K), 0.001, log = TRUE), tail(f$report$BBmsy, 1), f$report$r) }) }
Ks <- c(15e3, 20e3, 25e3, 30e3, 40e3, 50e3, 70e3, 90e3, 120e3)
for(ns in c(1, 4)){ pp <- prof(ns, Ks); d <- pp[1, ] - min(pp[1, ], na.rm = TRUE)
  cat("\nprofile ns =", ns, "\n"); print(round(rbind(K = Ks/1000, dNLL = d, BBmsy2021 = pp[2, ], r = pp[3, ]), 2)) }
