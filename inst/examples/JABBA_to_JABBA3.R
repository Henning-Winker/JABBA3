## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", fig.width = 9, fig.height = 5.5,
                      message = FALSE, warning = FALSE)


## ----packages-----------------------------------------------------------------
library(JABBA)        # JABBA 2 (JAGS)
library(JABBA3)
data(iccat, package = "JABBA3")
bet <- iccat$bet      # catch, cpue (2 longline indices) and se, 1950-2017


## ----jabba-fit----------------------------------------------------------------
jb <- build_jabba(catch = bet$catch, cpue = bet$cpue, se = bet$se, assessment = "BET", scenario = "Fox",
                  model.type = "Fox", sigma.est = FALSE, fixed.obsE = 0.01, psi.prior = c(1, 0.1), verbose = FALSE)
options(scipen = 0)                       # build_jabba() sets options(scipen = 999)
t.jabba <- system.time(jbf <- fit_jabba(jb, ni = 30000, nt = 5, nb = 10000, nc = 2, verbose = FALSE))[["elapsed"]]
jbf$pars[, 1:3]


## ----convert------------------------------------------------------------------
f.jabba <- jabba2jabba3(jbf, run = "JABBA")
f.jabba


## ----jabba-plots--------------------------------------------------------------
jbplot_trj(f.jabba)
jbplot_kobe(f.jabba)
jbplot_fits(f.jabba)
jbplot_residuals(f.jabba, type = "ordinary")
jbplot_runs(f.jabba, type = "ordinary")
jbplot_production(f.jabba)
jbplot_priorpost(f.jabba)


## ----refit--------------------------------------------------------------------
t.jabba3 <- system.time(f.jabba3 <- update(f.jabba, run = "JABBA3"))[["elapsed"]]
f.jabba3


## ----compare------------------------------------------------------------------
fits <- list(JABBA = f.jabba, JABBA3 = f.jabba3)
jbplot_trj(fits, type = c("B", "F", "B/Bmsy", "F/Fmsy"))
jb3_partable(fits, format = "wide")
data.frame(engine = c("JABBA (JAGS, 2 chains x 30,000)", "JABBA3 (RTMB)"), seconds = round(c(t.jabba, t.jabba3), 1))


## ----compare-diff-------------------------------------------------------------
d <- merge(f.jabba$ts[f.jabba$ts$quant %in% c("B/Bmsy", "F/Fmsy"), c("year", "quant", "est")],
           f.jabba3$ts[f.jabba3$ts$quant %in% c("B/Bmsy", "F/Fmsy"), c("year", "quant", "est")],
           by = c("year", "quant"), suffixes = c(".jabba", ".jabba3"))
aggregate(cbind(max.rel.diff = abs(est.jabba3/est.jabba - 1)) ~ quant, d, max)


## ----nuts, eval = requireNamespace("tmbstan", quietly = TRUE)-----------------
# t.nuts <- system.time(f.nuts <- update(f.jabba, run = "JABBA3 NUTS", engine = "nuts",
#                                        nuts = list(chains = 2, iter = 2000)))[["elapsed"]]
# jbplot_trj(list(JABBA = f.jabba, `JABBA3 NUTS` = f.nuts), type = c("B/Bmsy", "F/Fmsy"))
# round(t.nuts, 1)


## ----diffuse------------------------------------------------------------------
jb0 <- build_jabba(catch = bet$catch, cpue = bet$cpue, se = bet$se, assessment = "BET", scenario = "Fox igamma(0.001, 0.001)",
                   model.type = "Fox", sigma.est = FALSE, fixed.obsE = 0.01, psi.prior = c(1, 0.1),
                   igamma = c(0.001, 0.001), verbose = FALSE)
options(scipen = 0)
jbf0 <- fit_jabba(jb0, ni = 30000, nt = 5, nb = 10000, nc = 2, verbose = FALSE)
f0 <- jabba2jabba3(jbf0, run = "JABBA")
px <- list(`posterior x1` = update(jabba2jabba3(jbf0, sdb.inflate = 1, verbose = FALSE)),
           `posterior x1.5` = update(f0),
           `posterior x2` = update(jabba2jabba3(jbf0, sdb.inflate = 2, verbose = FALSE)),
           fixed = update(jabba2jabba3(jbf0, igamma.proxy = "fixed", verbose = FALSE)))
tab <- t(sapply(px, function(f){ m <- function(q){ a <- f$ts[f$ts$quant == q, "est"]; b <- f0$ts[f0$ts$quant == q, "est"]; max(abs(a/b - 1)) }
  c(sdb = f$pars["sdb", "est"], K = f$pars["K", "est"], `max rel. diff B/Bmsy` = m("B/Bmsy"), `max rel. diff F/Fmsy` = m("F/Fmsy")) }))
rbind(JABBA = c(f0$pars["sdb", "est"], f0$pars["K", "est"], 0, 0), tab)
jbplot_trj(c(list(JABBA = f0), px), type = c("B/Bmsy", "F/Fmsy"))


## ----list---------------------------------------------------------------------
jbfits <- list(Fox = jbf, `Fox igamma(0.001, 0.001)` = jbf0)
f.list <- lapply(jbfits, jabba2jabba3, verbose = FALSE)
jbplot_trj(f.list, type = c("B/Bmsy", "F/Fmsy"))
jb3_partable(f.list, format = "wide")

