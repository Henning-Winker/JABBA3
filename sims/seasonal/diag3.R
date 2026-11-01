suppressMessages(library(JABBA3, lib.loc = "NULL"))
src <- readLines("sim_seasonal.R"); i0 <- grep("^sim_om <-", src); i1 <- grep("^fitone <-", src) - 1
eval(parse(text = src[i0:i1]))
ft <- function(om, ns, it){ d <- jb3_data(om$catch, om$index, stock = "s", index.time = it)
  f <- suppressWarnings(jabba3(d, jb3_priors(d, r = c(1, 0.3)), jb3_model(BmsyK = 0.4, nseason = ns), osa = FALSE))
  bb <- f$ts$est[f$ts$quantity == "B/Bmsy"]; tail(log(bb) - log(om$truth$BBmsy), 1) }
for(it in c(0.5, 0)){
  r <- t(sapply(1:40, function(i){ om <- sim_om(i, w = rep(1, 4), fshare = rep(.25, 4), hmax = 1.5, itime = it)
    c(annual = ft(om, 1, it), seas_time_ok = ft(om, 4, it), seas_time_0 = ft(om, 4, 0)) }))
  cat("survey at", it, ":", paste(colnames(r), round(apply(r, 2, median), 3), collapse = "  "), "\n")
}
