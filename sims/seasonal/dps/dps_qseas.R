`%||%` <- function(a, b) if(is.null(a)) b else a
suppressMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); r <- e$ss3rep
cs <- aggregate(Obs ~ Yr + Seas, r$catch, sum); cs <- cs[order(cs$Yr, cs$Seas), ]
catch.s <- data.frame(year = cs$Yr, season = cs$Seas, catch = cs$Obs); catch.a <- aggregate(catch ~ year, catch.s, sum)
ci <- aggregate(cbind(Obs, SE) ~ Yr + Seas, r$cpue, mean); ci <- ci[order(ci$Yr, ci$Seas), ]
# survey split by season: one column per season (separate q), NA elsewhere
w <- reshape(ci[, c("Yr", "Seas", "Obs")], idvar = c("Yr", "Seas"), direction = "wide", timevar = "Seas")
mk <- function(v){ out <- data.frame(year = ci$Yr, season = ci$Seas)
  for(s in 2:4) out[[paste0("medits_s", s)]] <- ifelse(ci$Seas == s, ci[[v]], NA); out }
idx3 <- mk("Obs"); se3 <- mk("SE")
dS3 <- jb3_data(catch.s, idx3, se3, nseason = 4, stock = "DPS 12-16")
dA3 <- jb3_data(catch.a, idx3, se3, nseason = 4, stock = "DPS 12-16")
pr <- function(d){ p <- jb3_priors(d, r = c(1.07, 0.3), bkinit = c(0.9, 0.3)); p$sdb <- c(median = 0.15, log.sd = 0.5, on = 1); p }
runs <- list(q3_seascatch = list(d = dS3), q3_seascatch_prot = list(d = dS3, lag = 0.9, rho = 0.35),
             q3_anncatch = list(d = dA3), q3_anncatch_prot = list(d = dA3, lag = 0.9, rho = 0.35),
             q3_seascatch_sharedvar = list(d = dS3, sv = TRUE))
fits <- lapply(runs, function(x){
  md <- jb3_model(shape = "schaefer", nseason = 4, lag = x$lag %||% 0, rho = x$rho %||% 0, sdc = 0.2, sets.var = if(isTRUE(x$sv)) c(1, 1, 1) else NULL)
  tryCatch(suppressWarnings(jabba3(x$d, pr(x$d), md, osa = FALSE)), error = function(e){ message(conditionMessage(e)); NULL }) })
saveRDS(fits, "dps_fits_q3.rds")
tab <- do.call(rbind, lapply(names(fits), function(k){ f <- fits[[k]]; if(is.null(f)) return(data.frame(run = k, conv = "error"))
  g <- function(q) f$ts[f$ts$quantity == q, ]; bb <- g("B/Bmsy"); ff <- g("F/Fmsy"); n <- nrow(bb)
  data.frame(run = k, conv = f$convergence$tier, pdH = isTRUE(f$convergence$pdHess), r = round(f$report$r, 2), K = round(f$report$K*f$mucatch), sdb = round(f$report$sdb, 3),
             BBmsy2021 = sprintf("%.2f (%.2f-%.2f)", bb$est[n], bb$lci[n], bb$uci[n]), FFmsy2021 = sprintf("%.2f (%.2f-%.2f)", ff$est[n], ff$lci[n], ff$uci[n])) }))
print(tab, row.names = FALSE)
