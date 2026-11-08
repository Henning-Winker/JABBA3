# Thornback ray (Raja clavata) in Subarea 4 and divisions 3.a and 7.d (rjc.27.3a47d):
# SPiCT input and final assessment of the ICES benchmark WKBELASMO 2023 (ICES, 2023; doi:10.17895/ices.pub.22760042).
# Source: rjc_final.rdata (object `base`, the fitted SPiCT model of the final assessment), provided by H. Winker.
library(FLCore); library(FLRef)
e <- new.env(); load("rjc_final.rdata", envir = e); base <- e$base
inp0 <- base$inp
pr <- inp0$priors
on <- vapply(pr, function(x) if(is.list(x)) any(vapply(x, function(v) v[3] == 1, TRUE)) else length(x) > 2 && x[3] == 1, TRUE)
inp <- list(timeC = inp0$timeC, obsC = inp0$obsC, stdevfacC = inp0$stdevfacC,
            timeI = inp0$timeI, obsI = inp0$obsI, stdevfacI = inp0$stdevfacI,
            dteuler = inp0$dteuler, priors = pr[on])
spict.stk <- FLRef::spict2FLStockR(base, rel = TRUE)
rjc <- list(inp = inp, spict.stk = spict.stk,
            info = list(stock = "rjc.27.3a47d",
                        species = "Thornback ray (Raja clavata)",
                        area = "Subarea 4 and divisions 3.a and 7.d (North Sea, Skagerrak, Kattegat, eastern English Channel)",
                        source = "ICES (2023) WKBELASMO, final SPiCT assessment",
                        indices = c(Q14 = "NS-IBTS-Q1 and FR-CGFS-Q4 combined (biomass, >= 30 cm)",
                                    Q3 = "NS-IBTS-Q3, BTS-Eng-Q3 and BTS-BEL-Q3 combined (biomass, >= 30 cm)"),
                        lh = c(Linf = 123.62, k = 0.129, t0 = -1.154, L50 = 73.7, L95 = 82, a = 0.0045, b = 3.0961, M = 0.188, amax = 17)))
usethis::use_data(rjc, overwrite = TRUE)
