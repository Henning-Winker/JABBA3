fs <- commandArgs(TRUE)
lab <- c("trial_seasonal_asem_1_30.rds" = "none", "trial_seasonal_asem_1_30_protgt.rds" = "none", "trial_seasonal_asem_1_30_lagvar.rds" = "none",
         "trial_seasonal_asem_1_30_cap3.rds" = "0.30", "trial_seasonal_asem_1_30_cap25.rds" = "0.25")
x <- do.call(rbind, lapply(fs, function(f){ d <- readRDS(f); d$cap <- lab[[basename(f)]]; d }))
x <- x[x$year == 2025 & x$stock %in% c("pra.27.3a4a", "pil.27.8abd"), ]
x$ok <- x$tier %in% c("ok", "warning"); x$e <- log(x$est/x$tru_eb)
mods <- c("annual", "annual_protgt", "quarter", "quarter_lag50", "quarter_lagmax", "month_lag50", "month_lagmax")
x <- x[x$model %in% mods, ]
f <- function(d){ if(!nrow(d)) return(""); o <- d[d$ok, ]
  sprintf("%+.2f/%.2f/%.2f sdb %.2f%s", median(o$e), sqrt(mean(o$e^2)), mean(o$tru_eb >= o$lci & o$tru_eb <= o$uci, na.rm = TRUE),
          median(o$sdb), ifelse(mean(d$ok) < 1, sprintf(" (%.0f%%)", 100*mean(d$ok)), "")) }
for(st in c("pra.27.3a4a", "pil.27.8abd")){ cat("\n", st, " (median/RMSE/coverage, fitted sdb, convergence)\n")
  caps <- unique(x$cap)
  tab <- sapply(caps, function(cp) sapply(mods, function(m) f(x[x$stock == st & x$model == m & x$cap == cp, ])))
  print(noquote(tab)) }
