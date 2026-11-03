args <- commandArgs(TRUE); x <- do.call(rbind, lapply(args, readRDS))
ok <- x$tier %in% c("ok", "warning")
conv <- aggregate(ok ~ stock + model, data.frame(x[x$year == max(x$year), c("stock", "model")], ok = ok[x$year == max(x$year)]), mean)
x <- x[ok, ]
x$e_eb <- log(x$est/x$tru_eb); x$e_sb <- log(x$est/x$tru_sb)
x$phase <- cut(x$year, c(-Inf, 2004, 2014, Inf), c("down", "Fcut", "rebuild"))
mods <- c("annual", "annual_prot", "quarter_equal", "quarter", "quarter_prot", "month", "month_prot", "quarter_pulse", "quarter_pulse_prot", "month_pulse", "month_pulse_prot", "annual_protgt", "quarter_protgt", "month_protgt", "annual_ar", "quarter_ar", "month_ar", "quarter_lag50", "month_lag50", "quarter_lagmax", "month_lagmax")
term <- x[x$year == max(x$year), ]
s <- function(d) c(median = median(d$e_eb), rmse = sqrt(mean(d$e_eb^2)), x2 = mean(abs(d$e_eb) > log(2)),
                   cover = mean(d$tru_eb >= d$lci & d$tru_eb <= d$uci, na.rm = TRUE), vsSSB = median(d$e_sb))
tab <- do.call(rbind, lapply(split(term, list(term$stock, term$model), drop = TRUE), function(d) data.frame(stock = d$stock[1], model = d$model[1], t(s(d)))))
tab <- merge(tab, conv, by = c("stock", "model"))
ph <- aggregate(e_eb ~ stock + model + phase, x, median)
ph <- reshape(ph, idvar = c("stock", "model"), timevar = "phase", direction = "wide")
tab <- merge(tab, ph, by = c("stock", "model"))
tab$model <- factor(tab$model, mods); tab <- tab[order(tab$stock, tab$model), ]
names(tab) <- sub("e_eb.", "", names(tab))
cat("Terminal (2025) error vs true EB/EBmsy: median, RMSE, share off by > x2, 95% coverage; median vs SSB/SSBmsy; convergence;\n",
    "median error by phase (down 1981-2004, F cut 2005-14, rebuild 2015-25)\n", sep = "")
print(tab, digits = 2, row.names = FALSE)
cat("\nAll three stocks:\n")
all <- do.call(rbind, lapply(split(term, term$model), function(d) data.frame(model = d$model[1], t(s(d)))))
all$model <- factor(all$model, mods); print(all[order(all$model), ], digits = 2, row.names = FALSE)
