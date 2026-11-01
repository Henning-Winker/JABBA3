# Seasonal JABBA3 on the BioProcDyn age-structured OMs (annual OMs, recruitment at age 0 in January, index at mid-year)
# Annual vs seasonal (quarterly / monthly, equal-F split of the annual catch), each without (current practice) and with the
# protocol (production lag at amat + AR(1) rho 0.75). Priors as sims/R/jb3_treatments.R: r = OM asem r (log.sd 0.3),
# K prior off, bkinit = median OM EB depletion in 1981 (log.sd 0.3), FishLife sdb prior, Pella Bmsy/K 0.4, catch error.
# Usage: Rscript trial_seasonal_asem.R "1:30" cores [stocks]
suppressPackageStartupMessages({library(FLCore); library(FLBRP); library(parallel)})
args <- commandArgs(TRUE)
its <- eval(parse(text = if(length(args)) args[1] else "1:30")); cores <- if(length(args) > 1) as.integer(args[2]) else 2
st3 <- if(length(args) > 2) strsplit(args[3], ",")[[1]] else c("pra.27.3a4a", "her.27.25-2932", "pil.27.8abd")
jlib <- Sys.getenv("JB3LIB", "/home/claude/lib_new")
suppressPackageStartupMessages(library(JABBA3, lib.loc = jlib))
bpd <- Sys.getenv("BPD", "/home/claude/bioprocdyn")
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "flref_helpers.R"))

load(file.path(bpd, "data/sim9stks.rdata")); p8 <- readRDS(file.path(bpd, "output_procfloor/om_pil.27.8abd.rds"))
stks <- c(as.list(stksims), list("pil.27.8abd" = p8$stksim)); ix <- c(as.list(idxs), list("pil.27.8abd" = p8$idx))
bl <- c(as.list(brps), list("pil.27.8abd" = p8$brp)); rl <- c(as.list(rpts), list("pil.27.8abd" = p8$rpts))
year <- 1981

models <- list(
  annual         = list(ns = 1,  prot = FALSE),
  annual_prot    = list(ns = 1,  prot = TRUE),
  quarter_equal  = list(ns = 4,  prot = FALSE, alloc = "equal"),
  quarter        = list(ns = 4,  prot = FALSE),
  quarter_prot   = list(ns = 4,  prot = TRUE),
  month          = list(ns = 12, prot = FALSE),
  month_prot     = list(ns = 12, prot = TRUE),
  # recruitment on 1 January: half of the annual production in the step ending at the year boundary
  quarter_pulse      = list(ns = 4,  prot = FALSE, w = c(1, 1, 1, 3)),
  quarter_pulse_prot = list(ns = 4,  prot = TRUE,  w = c(1, 1, 1, 3)),
  month_pulse        = list(ns = 12, prot = FALSE, w = c(rep(1, 11), 11)),
  month_pulse_prot   = list(ns = 12, prot = TRUE,  w = c(rep(1, 11), 11)))
if(nzchar(Sys.getenv("MODELS"))) models <- models[strsplit(Sys.getenv("MODELS"), ",")[[1]]]

one <- function(it, st){
  s <- window(stks[[st]], start = year); idx <- window(ix[[st]], start = year); b <- bl[[st]]; rps <- rl[[st]]
  yrs <- as.numeric(dimnames(s)$year)
  itime <- mean(c(range(idx)[["startf"]], range(idx)[["endf"]]))          # survey timing within the year (0.5 here)
  inp <- data.frame(year = yrs, catch = c(iter(catch(s), it)), idx = c(iter(index(idx), it)))
  as <- asem2spm(b, quant = "vb"); ebmsy <- an(as["Bmsy"])
  b0 <- an(rps["B0"]); r0 <- an(rps["r"]); bk <- median(c(vb(stks[[st]])[, ac(year)]))/b0
  sbmsy <- an(refpts(b)["msy", "ssb"])
  tru_eb <- c(iter(vb(s), it))/ebmsy; tru_sb <- c(iter(ssb(s), it))/sbmsy
  amat <- as.numeric(dimnames(mat(b))$age)[which(c(mat(b)[, 1]) >= 0.5)[1]]
  d <- jb3_data(inp[c("year", "catch")], inp[c("year", "idx")], stock = st, index.time = itime)
  pr <- jb3_priors(d, r = c(r0, 0.3), K = c(b0, 1, 0), bkinit = c(bk, 0.3))
  do.call(rbind, lapply(names(models), function(nm){
    x <- models[[nm]]
    md <- jb3_model(BmsyK = 0.4, lag = if(x$prot) amat else 0, rho = if(x$prot) 0.75 else 0,
                    nseason = x$ns, catch.alloc = x$alloc %||% "F", season.w = x$w)
    t0 <- proc.time()[["elapsed"]]
    f <- tryCatch(suppressWarnings(jabba3(d, pr, md, osa = FALSE)), error = function(e) NULL)
    if(is.null(f)) return(data.frame(stock = st, iter = it, model = nm, year = yrs, est = NA, lci = NA, uci = NA,
                                     tru_eb = tru_eb, tru_sb = tru_sb, tier = "error", r = NA, sdb = NA, time = NA))
    g <- function(q, col = "est") f$ts[f$ts$quantity == q, col]
    data.frame(stock = st, iter = it, model = nm, year = yrs, est = g("B/Bmsy"), lci = g("B/Bmsy", "lci"), uci = g("B/Bmsy", "uci"),
               tru_eb = tru_eb, tru_sb = tru_sb, tier = f$convergence$tier, r = f$report$r, sdb = f$report$sdb,
               time = proc.time()[["elapsed"]] - t0)
  }))
}
`%||%` <- function(a, b) if(is.null(a)) b else a

grid <- expand.grid(st = st3, it = its, stringsAsFactors = FALSE)
t0 <- proc.time()[3]
res <- mclapply(seq_len(nrow(grid)), function(i) tryCatch(one(grid$it[i], grid$st[i]), error = function(e) conditionMessage(e)),
                mc.cores = cores)
bad <- !sapply(res, is.data.frame); if(any(bad)) print(unique(unlist(res[bad])))
out <- do.call(rbind, res[!bad])
cat("rows", nrow(out), "in", round(proc.time()[3] - t0), "s\n")
saveRDS(out, sprintf("trial_seasonal_asem_%s%s.rds", paste(range(its), collapse = "_"), Sys.getenv("TAG")))
