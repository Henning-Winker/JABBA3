# FLR version of the bigeye example input (JABBA iccat$bet), each series in its original length:
#   catch   FLQuant, 1950-2017 (total catch)
#   indices FLIndices of two FLIndexBiomass, each covering only the years with observations
#           (JR2_early 1959-1978, JR2_late 1979-2017); index.var left empty (NA) to show how to fill it
# Rscript data-raw/flr_betinp.R  (writes inst/extdata/flr_betinp.rds)
suppressMessages({library(FLCore); library(JABBA3)})
data(iccat); bet <- iccat$bet
x <- jbdat2flr(bet$catch, bet$cpue, bet$se)
crop <- function(ix){ v <- c(index(ix)); y <- as.numeric(dimnames(index(ix))$year)[is.finite(v)]
  ix <- window(ix, start = min(y), end = max(y)); index.var(ix)[] <- NA; ix }
flr_betinp <- list(catch = x$catch, indices = FLIndices(lapply(x$indices, crop)))
names(flr_betinp$indices) <- names(x$indices)
saveRDS(flr_betinp, "inst/extdata/flr_betinp.rds", version = 2)
str(lapply(flr_betinp$indices, function(i) range(as.numeric(dimnames(index(i))$year)))); print(range(as.numeric(dimnames(flr_betinp$catch)$year)))
