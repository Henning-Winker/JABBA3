lib <- commandArgs(TRUE)[1]
suppressMessages(library(JABBA3, lib.loc = lib))
data(iccat)
d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
pr <- jb3_priors(d, r = c(0.2, 0.5))
cfg <- list(fox = jb3_model(shape = "fox"), lag_ar = jb3_model(lag = 2, rho = 0.75), dyn = jb3_model(sdb.form = "dyn"),
            rw = jb3_model(rw = TRUE), bias_start = jb3_model(biascor = TRUE, start = "index", rho = 0.5),
            nocatcherr = jb3_model(catch.error = FALSE))
out <- lapply(cfg, function(m){ f <- suppressWarnings(jabba3(d, pr, m, osa = FALSE))
  c(nll = f$opt$objective, f$ts$est[f$ts$quantity %in% c("B/Bmsy","F/Fmsy","Catch","Process deviation")], f$refpts$est) })
saveRDS(out, commandArgs(TRUE)[2])
