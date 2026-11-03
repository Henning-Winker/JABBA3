suppressMessages(library(spict))
e <- new.env(); load("dps.12_16.SS3inp.rdata", envir = e); ss3rep <- e$ss3rep
df.catch <- ss3rep$catch
dfi <- aggregate(cbind(Obs, SE) ~ Yr + Time + Fleet_name + Fleet + Seas, ss3rep$cpue, mean); dfi <- dfi[order(dfi$Time), ]
dfc <- aggregate(cbind(Obs) ~ Time, df.catch, sum)
inp <- list(nseasons = 4, splineorder = 3)
inp$timeC <- floor(dfc$Time*4)/4; inp$obsC <- dfc$Obs
inp$stdevfacC <- rep(1, length(dfc$Obs)); inp$priors$logsdc <- c(log(0.2), 0.001, 1)
inp$timeI <- list(MeditsGrund = floor(dfi$Time*4)/4); inp$obsI <- list(MeditsGrund = dfi$Obs)
inp$stdevfacI <- list(dfi$SE/mean(dfi$SE, na.rm = TRUE))
oe <- mean(dfi$SE, na.rm = TRUE); inp$priors$logsdi <- list(MeditsGrund = c(log(oe), 0.3, 1))
r.pr <- c(1.072792, 0.3, 1); inp$priors$logr <- c(log(r.pr[1]) - 0.5*r.pr[2]^2, r.pr[2], r.pr[3])
inp$priors$logn <- c(log(2), 0.001, 1)
bk.pr <- c(0.9, 0.3, 1); inp$priors$logbkfrac <- c(log(bk.pr[1]) - bk.pr[2]^2/2, bk.pr[2], bk.pr[3])
pe <- c(0.15, 0.5, 1); inp$priors$logsdb <- c(log(pe[1]) - 0.5*pe[2]^2, pe[2], pe[3])
fdevs <- c(4, 0.5, 1); inp$priors$logsdf <- c(log(fdevs[1]) - 0.5*fdevs[2]^2, fdevs[[2]], fdevs[3])
inp$dteuler <- 1/4
ref <- fit.spict(inp)
cat("convergence", ref$opt$convergence, " pdHess", ref$pdHess, "\n")
g <- function(nm) get.par(nm, ref, exp = TRUE)
cat(sprintf("K %.0f  r %.2f  sdb %.3f  sdf %.2f  sdi %.3f  q %.2e\n", g("logK")[2], g("logr")[2], g("logsdb")[2], g("logsdf")[2], g("logsdi")[2], g("logq")[2]))
bb <- get.par("logBBmsy", ref, exp = TRUE); ff <- get.par("logFFmsy", ref, exp = TRUE); tt <- ref$inp$time
for(y in c(2019, 2020, 2021, 2022)){ i <- which.min(abs(tt - y)); cat(sprintf("t=%.2f  B/Bmsy %.2f (%.2f-%.2f)  F/Fmsy %.2f (%.2f-%.2f)\n", tt[i], bb[i, 2], bb[i, 1], bb[i, 3], ff[i, 2], ff[i, 1], ff[i, 3])) }
saveRDS(list(time = tt, bb = bb, ff = ff, K = g("logK")[2], r = g("logr")[2]), "spict_ref.rds")
# predicted vs observed catch (annual sums)
cp <- get.par("logCpred", ref, exp = TRUE)[, 2]; tc <- ref$inp$timeCpred
oc <- ref$inp$obsC; to <- ref$inp$timeC
yr <- floor(to + 1e-6); yp <- floor(tc + 1e-6)
ao <- tapply(oc, yr, sum); ap <- tapply(cp, yp, sum)
k <- intersect(names(ao), names(ap))
d <- data.frame(year = as.numeric(k), obs = as.numeric(ao[k]), pred = as.numeric(ap[k])); d$ratio <- d$pred/d$obs
print(d[d$year %in% c(1970, 1980, 1990, 2000, 2005, 2010, 2015, 2018, 2019, 2020, 2021), ], digits = 3, row.names = FALSE)
cat("median pred/obs 1969-2021:", round(median(d$ratio), 3), "  2010-2021:", round(median(d$ratio[d$year >= 2010]), 3), "\n")
