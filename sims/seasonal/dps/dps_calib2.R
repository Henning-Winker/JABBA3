`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
src <- readLines("dps_calib.R"); eval(parse(text = src[3:13]))
for(err in c("fixedSE", "spict")) for(cap in c(NA, 0.3, 0.22)){
  md <- jb3_model(shape = "schaefer", nseason = 4, sdc = 0.2, obs.error = if(err == "spict") "spict" else "additive", sdi.est = err != "fixedSE",
                  sdb.max = if(is.na(cap)) NULL else cap)
  f <- suppressWarnings(jabba3(dS, mk(err), md, osa = FALSE, nstart = 6))
  bb <- f$ts[f$ts$quantity == "B/Bmsy", ]; ff <- f$ts[f$ts$quantity == "F/Fmsy", ]; n <- nrow(bb)
  cat(sprintf("quarterly %-7s cap %-4s: %-4s K %6.0f r %.2f sdb %.3f sdi %.3f | B/Bmsy 2021 %.2f (%.2f-%.2f) F/Fmsy %.2f (%.2f-%.2f)\n", err, cap, f$convergence$tier,
      f$report$K*f$mucatch, f$report$r, f$report$sdb, f$report$sdi[1], bb$est[n], bb$lci[n], bb$uci[n], ff$est[n], ff$lci[n], ff$uci[n])) }
