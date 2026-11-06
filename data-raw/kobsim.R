# KOBsim: simulated silver kob (Argyrosomus inodorus) data set from the JABBA-Select example
# (github.com/Henning-Winker/JABBA-SELECTbeta, KOBsim_example; Winker et al. 2020, Fisheries Research 222: 105355).
# Age-structured operating model adapted from CCSRA (J. Thorson): 40 years, sigmaR = 0.6, minimum size limit raised in
# year 26 (length at 50% selectivity 300 -> 500 mm). Two indices of exploitable biomass (in weight), split at the change
# in selectivity (A1.S1, A1.S2, A2.S1, A2.S2); catch split by selectivity period (S1, S2).
src <- "/home/claude/jabba-selectbeta/KOBsim_example"
rd <- function(f) read.csv(file.path(src, "KOBsim", f))
catch <- rd("catchKOBsim.csv"); cpue <- rd("cpueKOBsim.csv"); se <- rd("seKOBsim.csv"); selex <- rd("selexKOBsim.csv")
names(catch)[1] <- names(cpue)[1] <- names(se)[1] <- "year"
load(file.path(src, "DataList_Original.RData"))
om <- DataList_Original
Lof <- function(a) 1372*(1 - exp(-0.115*(a + 0.815)))
kobsim <- list(
  catch = catch, cpue = cpue, se = se,
  # true operating-model biology (mature from age 4: maturity logistic in length between the lengths at ages 3 and 4)
  lhpar = c(linf = 1372, k = 0.115, t0 = -0.815, a = 6e-6, b = 3.07, M = 0.18, s = 0.8,
            L50 = round((Lof(3) + Lof(4))/2, 1), L95 = round(Lof(4) - 0.1*(Lof(4) - Lof(3)), 1), amax = 20),
  selpar = list(S1 = c(SL50 = selex$S1[1], SL95 = selex$S1[2]), S2 = c(SL50 = selex$S2[1], SL95 = selex$S2[2])),
  sel.change = 26,
  truth = data.frame(year = catch$year, SB = om$SB_t, EB = om$Bexploit_t, F = om$F_t, catch = om$Cw_t),
  SB0 = om$SB0,
  em.priors = list(M = c(0.2, 0.25), h = c(0.8, 0.1)),   # JABBA-Select estimation model: M (mean, CV), h (mean, CV)
  units = "length mm, weight g; catch and indices as simulated")
usethis::use_data(kobsim, overwrite = TRUE)
