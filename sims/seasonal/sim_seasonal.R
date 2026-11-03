# Simulation self-test of the seasonal JABBA3 prototype
# OM: seasonal Pella-Tomlinson SPM (4 seasons), production pulse in Q2, fishing concentrated in Q2-Q3,
# AR(1) process deviations, survey at the start of Q3 (time 0.5), catch reported annually only.
# Fits: annual (JABBA standard), seasonal equal-catch split, seasonal equal-F split, seasonal equal-F + known pulse.
args <- commandArgs(TRUE)
lib <- if(length(args)) args[1] else NULL
suppressMessages(library(JABBA3, lib.loc = lib))
`%||%` <- function(a, b) if(is.null(a)) b else a
niter <- if(length(args) > 1) as.integer(args[2]) else 20

sim_om <- function(seed, ns = 4, years = 1981:2025, r = 1.0, K = 1000, BmsyK = 0.4, w = c(0, 3, 1, 0),
                   fshare = c(0.1, 0.4, 0.4, 0.1), sdb = 0.2, rho = 0.5, sdi = 0.2, itime = 0.5, hmax = 2.5){
  set.seed(seed)
  m <- uniroot(function(m) m^(-1/(m - 1)) - BmsyK, c(1.01, 10))$root
  Hmsy <- r/(m - 1)*(1 - 1/m)
  Ny <- length(years); N <- Ny*ns; w <- w/mean(w)
  # annual harvest relative to Hmsy: fish-down to 2.5 x, then reduce to 0.7 x
  Hrel <- approx(c(1, 25, 30, Ny), c(0.3, hmax, hmax, 0.7), seq_len(Ny))$y
  rhos <- rho^(1/ns); f <- ns; for(k in 1:(ns - 1)) f <- f + 2*(ns - k)*rhos^k
  sds <- sdb/sqrt(f)
  P <- numeric(N + 1); P[1] <- 0.9; Cs <- numeric(N); d <- 0
  for(t in seq_len(N)){
    y <- (t - 1) %/% ns + 1; s <- (t - 1) %% ns + 1
    h <- min(Hrel[y]*Hmsy*fshare[s], 0.9)           # harvest fraction of the step
    Cs[t] <- h*P[t]*K
    pm <- max(P[t] + w[s]/ns*r/(m - 1)*P[t]*(1 - P[t]^(m - 1)) - Cs[t]/K, 0.001)
    d <- rhos*d + sqrt(1 - rhos^2)*rnorm(1, 0, sds)
    P[t + 1] <- pm*exp(d)
  }
  P <- P[1:N]
  yr1 <- (seq_len(Ny) - 1)*ns + 1
  Cy <- tapply(Cs, rep(seq_len(Ny), each = ns), sum)
  st <- yr1 + round(itime*ns)
  I <- 0.001*P[st]*K*exp(rnorm(Ny, 0, sdi))
  list(catch = data.frame(year = years, catch = as.numeric(Cy)), index = data.frame(year = years, survey = I),
       truth = data.frame(year = years, BBmsy = P[yr1]/BmsyK, BBmsy_mid = P[st]/BmsyK), m = m)
}

fitone <- function(om, ns, alloc = "F", w = NULL){
  d <- jb3_data(om$catch, om$index, stock = "sim", index.time = 0.5)
  pr <- jb3_priors(d, r = c(1.0, 0.3))
  md <- jb3_model(BmsyK = 0.4, nseason = ns, catch.alloc = alloc, season.w = w)
  f <- tryCatch(suppressWarnings(jabba3(d, pr, md, osa = FALSE)), error = function(e) NULL)
  if(is.null(f)) return(rep(NA, 3))
  bb <- f$ts[f$ts$quantity == "B/Bmsy", ]
  e <- log(bb$est) - log(om$truth$BBmsy)
  inci <- om$truth$BBmsy[nrow(bb)] >= bb$lci[nrow(bb)] & om$truth$BBmsy[nrow(bb)] <= bb$uci[nrow(bb)]
  c(terminal = tail(e, 1), rmse = sqrt(mean(e^2)), cover = inci)
}

runs <- list(annual = list(ns = 1), equal = list(ns = 4, alloc = "equal"), F = list(ns = 4, alloc = "F"),
             F_pulse = list(ns = 4, alloc = "F", w = c(0, 3, 1, 0)))
res <- do.call(rbind, lapply(seq_len(niter), function(i){
  om <- sim_om(i)
  do.call(rbind, lapply(names(runs), function(k){ x <- runs[[k]]
    v <- fitone(om, x$ns, x$alloc %||% "F", x$w); data.frame(iter = i, run = k, terminal = v[1], rmse = v[2], cover = v[3]) }))
}))
saveRDS(res, "sim_seasonal_res.rds")
agg <- do.call(rbind, lapply(split(res, res$run), function(x) data.frame(run = x$run[1], conv = mean(!is.na(x$terminal)),
  terminal_median = median(x$terminal, na.rm = TRUE), rmse_median = median(x$rmse, na.rm = TRUE), coverage = mean(x$cover, na.rm = TRUE))))
print(agg[match(names(runs), agg$run), ], digits = 3, row.names = FALSE)
