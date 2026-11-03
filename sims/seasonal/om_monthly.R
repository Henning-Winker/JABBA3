# Monthly operating models (close to continuous dynamics) for the seasonal JABBA3 prototype
# Annual catch only; survey at mid-year (start of July); fishing concentrated in months 4-9 (80% of annual F).
# Fits: annual, quarterly equal-catch, quarterly equal-F, quarterly equal-F with production weights matching the OM, monthly equal-F.
suppressMessages(library(JABBA3, lib.loc = "NULL"))
src <- readLines("sims/seasonal/sim_seasonal.R"); i0 <- grep("^sim_om <-", src); i1 <- grep("^fitone <-", src) - 1
eval(parse(text = src[i0:i1]))
niter <- 30
ft <- function(om, ns, alloc = "F", w = NULL){
  d <- jb3_data(om$catch, om$index, stock = "s", index.time = 0.5)
  f <- tryCatch(suppressWarnings(jabba3(d, jb3_priors(d, r = c(om$r, 0.3)),
         jb3_model(BmsyK = 0.4, nseason = ns, catch.alloc = alloc, season.w = w), osa = FALSE)), error = function(e) NULL)
  if(is.null(f) || f$convergence$tier == "fail") return(c(NA, NA, NA))
  bb <- f$ts[f$ts$quantity == "B/Bmsy", ]; e <- log(bb$est) - log(om$truth$BBmsy); n <- nrow(bb)
  c(tail(e, 1), sqrt(mean(e^2)), om$truth$BBmsy[n] >= bb$lci[n] & om$truth$BBmsy[n] <= bb$uci[n])
}
fs12 <- c(rep(0.2/6, 3), rep(0.8/6, 6), rep(0.2/6, 3))
prod <- list(continuous = list(w12 = rep(1, 12), w4 = NULL),
             one_pulse  = list(w12 = c(0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0), w4 = c(0, 1, 0, 0)),     # recruitment Apr-Jun (Q2)
             two_pulses = list(w12 = c(0, 0, 0, 2, 2, 0, 0, 0, 0, 2, 2, 0), w4 = c(0, 1, 0, 1)))     # Apr-May (Q2) and Oct-Nov (Q4)
out <- list()
for(pn in names(prod)) for(r in c(0.3, 1, 2)){
  res <- lapply(seq_len(niter), function(i){
    om <- sim_om(i, ns = 12, r = r, w = prod[[pn]]$w12, fshare = fs12, hmax = 1.5, itime = 0.5); om$r <- r
    fits <- list(annual = ft(om, 1), q_equal = ft(om, 4, "equal"), q_F = ft(om, 4, "F"),
                 q_F_w = if(is.null(prod[[pn]]$w4)) rep(NA, 3) else ft(om, 4, "F", prod[[pn]]$w4), m_F = ft(om, 12, "F"))
    data.frame(prod = pn, r = r, iter = i, model = names(fits), terminal = sapply(fits, `[`, 1), rmse = sapply(fits, `[`, 2),
               cover = sapply(fits, `[`, 3), minBB = min(om$truth$BBmsy))
  })
  out[[paste(pn, r)]] <- do.call(rbind, res)
  cat(pn, r, "done\n")
}
res <- do.call(rbind, out); saveRDS(res, "sims/seasonal/om_monthly_res.rds")
agg <- aggregate(cbind(terminal, rmse, cover) ~ prod + r + model, res, function(x) median(x), na.action = na.pass)
agg$cover <- aggregate(cover ~ prod + r + model, res, mean, na.action = na.omit)$cover
agg$conv <- aggregate(terminal ~ prod + r + model, res, function(x) mean(!is.na(x)), na.action = na.pass)$terminal
agg$model <- factor(agg$model, c("annual", "q_equal", "q_F", "q_F_w", "m_F"))
agg <- agg[order(agg$prod, agg$r, agg$model), ]
print(agg, digits = 2, row.names = FALSE)
cat("min true B/Bmsy (median over iterations):\n"); print(aggregate(minBB ~ prod + r, res, median), digits = 2)
