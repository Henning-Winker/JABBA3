# HYPOTHETICAL auxiliary series for bigeye tuna (iccat$bet$effort_hyp, iccat$bet$ffmsy_hyp).
# Not real data: simulated from the F trajectory of the JABBA3 bet example fit (Fox, r prior c(0.2, 0.5), both JPN LL
# indices) to illustrate auxiliary data in the user guide.
#   effort_hyp: E_t = F_t / (q0 x 1.02^(t - 1960)) x exp(eps), eps ~ N(0, 0.15^2), 1960-2017, q0 = 1e-4
#               (2% annual effort creep; units arbitrary, "relative effort")
#   ffmsy_hyp:  F/Fmsy_t x exp(eps), eps ~ N(0, 0.3^2), 1990-2009 only (a high-exploitation period), se = 0.3
library(JABBA3)
load("data/iccat.rda")
bet <- iccat$bet
d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
f <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"), osa = FALSE)
F <- f$ts[f$ts$quantity == "F", ]; FF <- f$ts[f$ts$quantity == "F/Fmsy", ]
set.seed(2025)
y <- 1960:2017; i <- match(y, F$year)
eff <- F$est[i]/(1e-4*1.02^(y - 1960))*exp(stats::rnorm(length(y), 0, 0.15))
iccat$bet$effort_hyp <- data.frame(year = y, effort = signif(eff, 4))
y2 <- 1990:2009; i2 <- match(y2, FF$year)
iccat$bet$ffmsy_hyp <- data.frame(year = y2, ffmsy = signif(FF$est[i2]*exp(stats::rnorm(length(y2), 0, 0.3)), 4), se = 0.3)
iccat$bet$info$hypothetical <- paste("effort_hyp and ffmsy_hyp are HYPOTHETICAL series simulated from the bet F trajectory",
  "(effort: 2% annual creep, lognormal noise SD 0.15; F/Fmsy: noise SD 0.3, 1990-2009); not real data")
save(iccat, file = "data/iccat.rda", compress = "xz")
