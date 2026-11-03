`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs); mc <- max(aggregate(catch ~ year, catch.s, sum)$catch)
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
dS <- jb3_data(catch.s, data.frame(year = ci$Yr, season = ci$Seas, medits = ci$Obs), data.frame(year = ci$Yr, season = ci$Seas, medits = ci$SE), nseason = 4)
oe <- mean(ci$SE)
mk <- function(err){ p <- jb3_priors(dS, r = c(1.072792*exp(-0.5*0.3^2), 0.3), bkinit = c(0.9*exp(-0.5*0.3^2), 0.3),
                                     K = c(4*mc, 0.833, 1)); p$sdb <- c(median = 0.15*exp(-0.5*0.5^2), log.sd = 0.5, on = 1)
  if(err == "spict") p$sdi <- c(median = oe, log.sd = 0.3, on = 1); p }
sp <- readRDS("spict_ref.rds"); cat(sprintf("K prior median 4 x max catch = %.0f\n", 4*mc))
res <- list()
for(err in c("additive", "fixedSE", "spict")) for(ns in c(1, 4)) for(rho in c(0, 0.35)){
  md <- jb3_model(shape = "schaefer", nseason = ns, sdc = 0.2, rho = rho, obs.error = if(err == "spict") "spict" else "additive",
                  sdi.est = err != "fixedSE")
  f <- tryCatch(suppressWarnings(jabba3(dS, mk(err), md, osa = FALSE, nstart = 6)), error = function(e) NULL)
  if(is.null(f)){ cat(err, ns, rho, "error\n"); next }
  bb <- f$ts[f$ts$quantity == "B/Bmsy", ]; ff <- f$ts[f$ts$quantity == "F/Fmsy", ]; n <- nrow(bb)
  res[[length(res) + 1]] <- data.frame(index_error = err, model = if(ns == 1) "annual" else "quarterly", rho = rho, conv = f$convergence$tier,
    K = round(f$report$K*f$mucatch), r = round(f$report$r, 2), sdb = round(f$report$sdb, 3), sdi = round(f$report$sdi[1], 3),
    BBmsy2021 = sprintf("%.2f (%.2f-%.2f)", bb$est[n], bb$lci[n], bb$uci[n]), FFmsy2021 = sprintf("%.2f (%.2f-%.2f)", ff$est[n], ff$lci[n], ff$uci[n]))
}
tab <- do.call(rbind, res); print(tab, row.names = FALSE)
cat(sprintf("\nSPiCT: K %.0f r %.2f | B/Bmsy 2021 0.76 (0.43-1.33) F/Fmsy 1.12 (0.62-2.03) [Jan 1; 0.87 mid-year]\nSS3: SSB/SSBmsy 0.70, F/Fmsy 0.78\n", sp$K, sp$r))
saveRDS(tab, "dps_calib_tab.rds")
