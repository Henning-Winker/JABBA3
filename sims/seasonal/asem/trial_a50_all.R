# 50%-maturity lag check on all nine case-study OMs (annual integer lag vs annual lag at round(a50) vs seasonal fractional a50 lag)
suppressPackageStartupMessages({library(FLCore); library(FLBRP); library(parallel)})
args <- commandArgs(TRUE); its <- eval(parse(text = args[1])); cores <- as.integer(args[2])
suppressPackageStartupMessages(library(JABBA3, lib.loc = "/home/claude/lib_new"))
bpd <- "/home/claude/bioprocdyn"; source("/home/claude/work/seasonal_om/flref_helpers.R")
`%||%` <- function(a, b) if(is.null(a)) b else a
load(file.path(bpd, "data/sim9stks.rdata"))
p8 <- readRDS(file.path(bpd, "output_procfloor/om_pil.27.8abd.rds")); oa <- readRDS(file.path(bpd, "output_procfloor/om_anf.27.3a46.rds"))
stks <- c(as.list(stksims), list("pil.27.8abd" = p8$stksim, "anf.27.3a46" = oa$stksim)); ix <- c(as.list(idxs), list("pil.27.8abd" = p8$idx, "anf.27.3a46" = oa$idx))
bl <- c(as.list(brps), list("pil.27.8abd" = p8$brp, "anf.27.3a46" = oa$brp)); rl <- c(as.list(rpts), list("pil.27.8abd" = p8$rpts, "anf.27.3a46" = oa$rpts))
st9 <- c("pil.27.8abd", "pra.27.3a4a", "her.27.25-2932", "hke.27.8c9a", "ple.27.420", "hom.27.2a4a5b6a7a-ce-k8", "bss.27.8ab", "cod.27.5a", "anf.27.3a46")
year <- 1981
# per-stock biology: integer amat (first age mat >= 0.5), interpolated a50 maturity, generation time and table rho
bio <- do.call(rbind, lapply(st9, function(st){
  b <- bl[[st]]; a <- as.numeric(dimnames(mat(b))$age); m <- c(mat(b)[, 1])
  k <- which(m >= 0.5)[1]; a50 <- if(k == 1) a[1] else a[k - 1] + (0.5 - m[k - 1])/(m[k] - m[k - 1])
  gt <- tryCatch(jb3_leslie(iter(stks[[st]], 1), h = an(rl[[st]]["s"]))[["gt"]], error = function(e) NA)
  rho <- if(is.na(gt)) 0.75 else if(gt <= 3) 0.35 else if(gt <= 6) 0.6 else if(gt <= 10) 0.75 else 0.8
  data.frame(stock = st, amat = a[k], a50 = round(a50, 2), a50y = max(round(a50), 1), gt = round(gt, 1), rho = rho)
}))
print(bio, row.names = FALSE)
models <- list(annual = list(ns = 1, lag = "none", rho = FALSE), annual_prot = list(ns = 1, lag = "amat", rho = TRUE),
               annual_lag50 = list(ns = 1, lag = "a50y", rho = TRUE), quarter_lag50 = list(ns = 4, lag = "a50", rho = TRUE),
               month_lag50 = list(ns = 12, lag = "a50", rho = TRUE))
one <- function(it, st){
  s <- window(stks[[st]], start = year); idx <- window(ix[[st]], start = year); b <- bl[[st]]; rps <- rl[[st]]; bi <- bio[bio$stock == st, ]
  yrs <- as.numeric(dimnames(s)$year); itime <- mean(c(range(idx)[["startf"]], range(idx)[["endf"]]))
  inp <- data.frame(year = yrs, catch = c(iter(catch(s), it)), idx = c(iter(index(idx), it)))
  ebmsy <- an(asem2spm(b, quant = "vb")["Bmsy"]); b0 <- an(rps["B0"]); r0 <- an(rps["r"])
  bk <- median(c(vb(stks[[st]])[, ac(year)]))/b0; sbmsy <- an(refpts(b)["msy", "ssb"])
  tru_eb <- c(iter(vb(s), it))/ebmsy; tru_sb <- c(iter(ssb(s), it))/sbmsy
  d <- jb3_data(inp[c("year", "catch")], inp[c("year", "idx")], stock = st, index.time = itime)
  pr <- jb3_priors(d, r = c(r0, 0.3), K = c(b0, 1, 0), bkinit = c(bk, 0.3))
  do.call(rbind, lapply(names(models), function(nm){ x <- models[[nm]]
    lg <- switch(x$lag, none = 0, amat = bi$amat, a50y = bi$a50y, a50 = bi$a50)
    md <- jb3_model(BmsyK = 0.4, lag = lg, rho = if(x$rho) bi$rho else 0, nseason = x$ns)
    f <- tryCatch(suppressWarnings(jabba3(d, pr, md, osa = FALSE)), error = function(e) NULL)
    if(is.null(f)) return(data.frame(stock = st, iter = it, model = nm, year = yrs, est = NA, lci = NA, uci = NA, tru_eb = tru_eb, tru_sb = tru_sb, tier = "error", r = NA, sdb = NA))
    g <- function(q, col = "est") f$ts[f$ts$quantity == q, col]
    data.frame(stock = st, iter = it, model = nm, year = yrs, est = g("B/Bmsy"), lci = g("B/Bmsy", "lci"), uci = g("B/Bmsy", "uci"),
               tru_eb = tru_eb, tru_sb = tru_sb, tier = f$convergence$tier, r = f$report$r, sdb = f$report$sdb) }))
}
grid <- expand.grid(st = st9, it = its, stringsAsFactors = FALSE)
t0 <- proc.time()[3]
res <- mclapply(seq_len(nrow(grid)), function(i) tryCatch(one(grid$it[i], grid$st[i]), error = function(e) conditionMessage(e)), mc.cores = cores)
bad <- !sapply(res, is.data.frame); if(any(bad)) print(unique(unlist(res[bad])))
out <- do.call(rbind, res[!bad]); cat("rows", nrow(out), "in", round(proc.time()[3] - t0), "s\n")
saveRDS(list(res = out, bio = bio), sprintf("trial_a50_all_%s.rds", paste(range(its), collapse = "_")))
