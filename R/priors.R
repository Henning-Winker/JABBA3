#' Priors for JABBA3
#'
#' One named list of priors, each \code{c(median, log.sd, on)}: lognormal with the given median and log-scale sd, and an
#' on/off switch (\code{on = 0} gives a flat prior on the log scale). The defaults follow JABBA
#' (\code{build_jabba()}), except the process error prior.
#'
#' @param data optional \code{jb3data}, used for the default K prior (8 x max catch, CV 1, as JABBA)
#' @param r intrinsic growth rate, on by default
#' @param K carrying capacity prior c(median, log.sd[, on]); \code{NULL} or \code{"jabba"} = 8 x max catch with log.sd 0.833
#'   (CV 1; the JABBA default, \code{build_jabba(K.prior = NULL)}); \code{"jabbaN"} = N x max catch, CV 1 (e.g. \code{"jabba4"})
#'   \code{"half6"} (or \code{"half"}) = one-sided empirical prior: flat below 6 x max catch and a half-normal on log K above
#'   it with log-SD 0.55 (median 5.9 and between-stock log-SD 0.55 of K / max catch in 73 ICES stocks; the data supply the
#'   lower bound, so only the upper side needs a prior); \code{"halfN"} = the same at N x max catch (e.g. \code{"half4"} for
#'   short generation times or series that mostly reflect depleted or rebuilding states); \code{"half+r"} = at the median
#'   predicted from r, \code{jb3_kr(r)} (3.04 x r^-0.583 x max catch, r = median of the r prior). Implemented through \code{K.max} (the K prior is off). The
#'   location is a median: for an expected value E (e.g. JABBA's 8 x with CV 1), use E x exp(-log.sd^2/2) (8 x gives 5.7 x)
#' @param bkinit initial depletion B1/K (JABBA: \code{psi}); default c(0.9, 0.246) (JABBA: CV 0.25)
#' @param sdb process error SD: \code{"fishlife"} (default) = LN(0.80 x sigR x r^0.45, 0.3) with r the median of the
#'   r prior (FishLife relation for the SD at Bmsy); or \code{c(median, log.sd, on)}; or \code{"igamma"} for JABBA's
#'   diffuse inverse-gamma(0.001, 0.001) on the variance (not recommended for the posterior mode).
#' @param sigR recruitment SD used by \code{sdb = "fishlife"} (0.5 when unknown)
#' @param q catchability: \code{NULL} = flat on log q for every q (default); otherwise a matrix (or vector) with one row
#'   \code{c(median, log.sd, on)} per catchability group
#' @param m Pella-Tomlinson shape: off by default (m fixed from \code{BmsyK} in \code{\link{jb3_model}()});
#'   \code{c(median, log.sd, 1)} to estimate m
#' @param kdyn depletion dependence of the process error SD (\code{sdb.form = "dyn"}), \code{c(mean, sd, on)} normal
#'   on k, bounded to (-0.5, 2); \code{NULL} = life-history centre \code{max(0, 0.064 - 0.155 log r)}, sd 0.15
#' @param rho AR(1) coefficient of the process deviations, \code{c(mean, sd, on)} normal on rho in (0, 0.99); off by
#'   default (rho fixed by \code{jb3_model(rho = )})
#' @param sdi index observation SD prior (only with \code{obs.error = "spict"}), off by default: one prior
#'   \code{c(median, log.sd, on)} for all observation-variance groups, or a list (or 3-column matrix) with one prior per
#'   group (\code{sets.var} in \code{jb3_model()}; by default one group per index), as SPiCT's \code{logsdi} by index
#' @param status optional status priors: a data.frame(year, type, median, log.sd) with type one of "BK", "BBmsy",
#'   "FFmsy" (lower case and JABBA's "bk", "bbmsy", "ffmsy" accepted; year NA = terminal year of the fit), the output of
#'   \code{\link{jb3_bprior}()}, or a JABBA \code{b.prior} vector c(mean, cv, year, type) (the cv used as log-scale SD, as
#'   in JABBA; year NA = terminal year)
#' @param sdf SD of the random walk in log harvest rate (\code{jb3_model(frw = )}): lognormal prior c(median, log.sd, on),
#'   default LN(0.3, 1) (95\% range about 0.04-2). With \code{frw = "gaps"} the data carry little information on sdf (the
#'   prior sets it); with \code{"index"} or \code{"all"} and catch known, sdf is estimated from the harvest-rate path
#'   (typically 0.1-0.4). A prior centred far above that, e.g. SPiCT's c(4, 0.5), pulls the estimate up 1.3-2x.
#' @param creep normal prior c(mean, sd) on the annual effort creep, used with \code{jb3_model(effort.creep = "estimate")}
#' @param psi alias of \code{bkinit} (JABBA name)
#' @param bbmsyinit initial status as B/Bmsy instead of B/K: c(median, log.sd), lognormal on B/Bmsy in the first year (a
#'   status prior; the bkinit prior is switched off). With a fixed shape it is equivalent to bkinit = median x BmsyK, but it
#'   keeps the initial status the same across shapes, e.g. for a BmsyK sensitivity; needs \code{data}
#' @param sdc catch SD prior c(median, log.sd, on) when \code{jb3_model(sdc.est = TRUE)} (off by default: flat on log sdc)
#' @param maxC the maximum catch used for the K prior and the K.max cap: \code{"max"} (default; maximum annual total catch,
#'   as JABBA) or \code{"smooth"} (maximum of the catch smoothed with a normal kernel, bandwidth 3 years, as CMSY /
#'   \code{JABBA::cmsy2jabba()}; less sensitive to a single peak year)
#' @param K.max optional soft upper prior on K, to avoid implausibly large biomass when the data say little about
#'   the scale: no penalty below the cap, and a half-normal on log K above it. \code{"r"} sets the cap from r (the median of
#'   the r prior) at the 95th percentile of K / max catch across 72 ICES stocks (K = unfished exploitable biomass):
#'   cap = 7.6 x max catch x r^-0.587 (about 19 x max catch at r = 0.2, 12 x at 0.4, 7.6 x at 1); a number gives the cap
#'   as a multiple of the maximum catch; \code{c(multiple, log.sd)} also sets the SD of the half-normal (default 0.3).
#'   \code{NULL} (default) = no cap. Soft: 1.5 x the cap adds 0.9 to the objective, 2 x adds 2.7, 3 x adds 6.7. Not on by
#'   default because lightly exploited, low-productivity stocks (e.g. rays) can lie well above the r-based cap; instead,
#'   \code{jabba3()} warns when K exceeds it, and the \code{inits = "best"} search tries a start from a capped fit.
#' @return list of class \code{jb3priors}
#' @examples
#' jb3_priors(r = c(0.3, 0.3))
#' @export
jb3_priors <- function(data = NULL, r = c(0.2, 0.5, 1), K = NULL, bkinit = c(0.9, 0.246, 1), sdb = "fishlife", sigR = 0.5,
                       q = NULL, m = c(NA, 0.3, 0), kdyn = NULL, rho = c(0.7, 0.15, 0), sdi = c(0.2, 0.5, 0), status = NULL, psi = NULL, bbmsyinit = NULL,
                       creep = c(0.02, 0.02), sdf = c(0.3, 1, 1), K.max = NULL, sdc = c(0.1, 1, 0), maxC = c("max", "smooth")){
  maxC <- match.arg(maxC)
  args <- .fn_args()
  pr3 <- function(x, nm){
    x <- as.numeric(x); if(length(x) == 2) x <- c(x, 1)
    if(length(x) != 3) stop("prior ", nm, " must be c(median, log.sd[, on])", call. = FALSE)
    stats::setNames(x, c("median", "log.sd", "on"))
  }
  if(!is.null(psi)) bkinit <- psi
  r <- pr3(r, "r")
  mC <- if(!is.null(data)) .max_catch(data, maxC) else NA
  if(is.null(K)) K <- "jabba"
  if(is.character(K) && grepl("^half", K)){      # one-sided empirical K prior: flat below, half-normal on log K above
    loc <- if(identical(K, "half")) 6 else if(K %in% c("half+r", "half.r")) jb3_kr(r = as.numeric(r[1]))[["median"]] else
      if(grepl("^half[0-9.]+$", K)) as.numeric(sub("^half", "", K)) else
      stop("K: \"half\", \"halfN\" (e.g. \"half6\", \"half4\") or \"half+r\"", call. = FALSE)
    if(!is.null(K.max)) warning("K = \"", K, "\" replaces K.max", call. = FALSE)
    K.max <- c(loc, 0.55); K <- c(loc*mC, 0.55, 0)
  }
  if(is.character(K)){                           # "jabba" = 8 x max catch, CV 1 (JABBA default); "jabbaN" = N x max catch, CV 1
    mult <- if(grepl("^jabba[0-9.]+$", K)) as.numeric(sub("^jabba", "", K)) else if(identical(K, "jabba")) 8 else
      stop("K must be c(median, log.sd[, on]), \"jabba\" (8 x max catch, CV 1), \"jabbaN\" (N x max catch, CV 1), \"halfN\" or \"half+r\"", call. = FALSE)
    K <- c(mult*mC, cv2sd(1), 1)
  }
  K <- pr3(K, "K"); bkinit <- pr3(bkinit, "bkinit")
  sdb.type <- if(is.character(sdb)) match.arg(sdb, c("fishlife", "igamma")) else "lognormal"
  sdb <- if(sdb.type == "fishlife") pr3(c(0.80*sigR*r[["median"]]^0.45, 0.3, 1), "sdb") else
    if(sdb.type == "igamma") c(median = NA, log.sd = NA, on = 1) else pr3(sdb, "sdb")
  if(!is.null(q)){ q <- matrix(q, ncol = 3, dimnames = list(NULL, c("median", "log.sd", "on"))) }
  m <- pr3(m, "m")
  if(is.null(kdyn)) kdyn <- c(if(r[["median"]] > 1.5) 0 else max(0, 0.064 - 0.155*log(r[["median"]])), 0.15, 1)
  kdyn <- stats::setNames(as.numeric(c(kdyn, 1)[1:3]), c("mean", "sd", "on"))
  sdi <- if(is.list(sdi) || is.matrix(sdi)) .pr_rows(sdi, "sdi") else pr3(sdi, "sdi"); sdf <- pr3(sdf, "sdf")
  rho <- stats::setNames(as.numeric(c(rho, 1)[1:3]), c("mean", "sd", "on"))
  status <- .status_df(status)
  if(!is.null(bbmsyinit)){                       # initial status as B/Bmsy instead of B/K: status prior in the first year, bkinit off
    if(is.null(data)) stop("bbmsyinit needs data (to set the first year)", call. = FALSE)
    b <- as.numeric(bbmsyinit); if(length(b) < 2) stop("bbmsyinit must be c(median, log.sd)", call. = FALSE)
    status <- rbind(status, data.frame(year = min(data$years), type = "BBmsy", median = b[1], log.sd = b[2]))
    bkinit[["on"]] <- 0
  }
  K.max <- if(is.null(K.max)) c(mult = NA, log.sd = NA, on = 0) else if(identical(K.max, "r"))
    c(mult = 7.6*r[["median"]]^-0.587, log.sd = 0.3, on = 1) else {
      if(!is.numeric(K.max) || K.max[1] <= 0) stop("K.max must be \"r\", a multiple of the maximum catch, or c(multiple, log.sd)", call. = FALSE)
      c(mult = K.max[1], log.sd = if(length(K.max) > 1) K.max[2] else 0.3, on = 1) }
  structure(list(r = r, K = K, bkinit = bkinit, sdb = sdb, sdb.type = sdb.type, sigR = sigR, q = q, m = m, kdyn = kdyn, rho = rho,
                 sdi = sdi, sdf = sdf, sdc = stats::setNames(as.numeric(sdc), c("median", "log.sd", "on")), status = status, creep = stats::setNames(as.numeric(creep[1:2]), c("mean", "sd")), K.max = K.max,
                 maxC = maxC, maxC.value = mC),
            class = "jb3priors", args = args)
}

#' Status priors (JABBA b.prior)
#'
#' Prior on stock status in a given year: B/Bmsy, B/K or F/Fmsy, lognormal with the given median and CV (as JABBA's
#' \code{b.prior = c(mean, cv, year, type)}). Several priors are combined with \code{rbind()} and passed to
#' \code{jb3_priors(status = )}.
#' @param median median of the prior (e.g. 0.5 for B/Bmsy)
#' @param cv coefficient of variation (converted to a log-scale SD, sqrt(log(1 + cv^2)))
#' @param year year(s) of the prior; NULL or NA = terminal year of the fit
#' @param type "bbmsy" (default), "bk" or "ffmsy"
#' @param log.sd log-scale SD, instead of \code{cv}
#' @return data.frame(year, type, median, log.sd)
#' @examples
#' jb3_bprior(0.7, 0.2, type = "bbmsy")                # terminal-year B/Bmsy prior
#' rbind(jb3_bprior(0.5, 0.2, 2010), jb3_bprior(1.2, 0.3, 2010, "ffmsy"))
#' @export
jb3_bprior <- function(median, cv = 0.2, year = NULL, type = c("bbmsy", "bk", "ffmsy"), log.sd = NULL){
  type <- match.arg(type)
  sd <- log.sd %||% cv2sd(cv)
  data.frame(year = if(is.null(year)) NA_real_ else as.numeric(year), type = c(bbmsy = "BBmsy", bk = "BK", ffmsy = "FFmsy")[[type]],
             median = as.numeric(median), log.sd = as.numeric(sd))
}

.status_df <- function(status){
  if(is.null(status)) return(NULL)
  if(is.atomic(status) && !is.matrix(status)){                      # JABBA b.prior c(mean, cv, year, type)
    if(length(status) < 4) stop("b.prior vector must be c(mean, cv, year, type)", call. = FALSE)
    yr <- suppressWarnings(as.numeric(status[3]))
    status <- jb3_bprior(as.numeric(status[1]), year = if(is.na(yr)) NULL else yr, type = tolower(status[4]),
                         log.sd = as.numeric(status[2]))              # JABBA uses the cv as the log-scale SD
  }
  status <- as.data.frame(status)
  if(!all(c("year", "type", "median", "log.sd") %in% names(status))) stop("status priors need year, type, median, log.sd", call. = FALSE)
  map <- c(bk = "BK", bbmsy = "BBmsy", ffmsy = "FFmsy")
  tp <- as.character(status$type); lw <- tolower(tp); tp[lw %in% names(map)] <- map[lw[lw %in% names(map)]]
  status$type <- tp
  if(!all(status$type %in% c("BK", "BBmsy", "FFmsy"))) stop("status prior type must be BK, BBmsy or FFmsy", call. = FALSE)
  status$year <- as.numeric(status$year)
  status
}

#' @export
print.jb3priors <- function(x, ...){
  sdi <- if(is.matrix(x$sdi)){ s <- x$sdi; rownames(s) <- paste0("sdi.", rownames(s) %||% seq_len(nrow(s))); s } else rbind(sdi = x$sdi)
  tab <- rbind(r = x$r, K = x$K, bkinit = x$bkinit, sdb = x$sdb, m = x$m, sdi)
  cat("JABBA3 priors (lognormal: median, log.sd; on = 1)\n")
  print(signif(tab, 3))
  cat("sdb:", switch(x$sdb.type, fishlife = sprintf("FishLife centre 0.80 x sigR (%.2f) x r^0.45", x$sigR),
                     igamma = "diffuse inverse-gamma(0.001, 0.001) on the variance (JABBA legacy)", lognormal = "user lognormal"), "\n")
  cat("q:", if(is.null(x$q)) "flat on log q" else paste(nrow(x$q), "lognormal priors"), "\n")
  cat(sprintf("kdyn (dynamic sdb only): normal(%.2f, %.2f)%s\n", x$kdyn[1], x$kdyn[2], if(x$kdyn[3] == 0) " fixed at mean" else ""))
  if(!is.null(x$maxC.value) && is.finite(x$maxC.value)) cat(sprintf("max catch (%s): %s\n", x$maxC %||% "max", format(signif(x$maxC.value, 4), big.mark = ",")))
  if(!is.null(x$K.max) && x$K.max[3] == 1) cat(sprintf("K.max: soft upper prior at %.1f x max catch (half-normal on log K above, sd %.2f)\n", x$K.max[1], x$K.max[2]))
  if(!is.null(x$status)){ cat("status priors:\n"); print(x$status) }
  invisible(x)
}

# one prior c(median, log.sd, on) per row, from a list of vectors or a 3-column matrix (rows named after list names)
.pr_rows <- function(x, nm){
  if(is.list(x)) x <- do.call(rbind, lapply(x, function(v){ v <- as.numeric(v); if(length(v) == 2) v <- c(v, 1)
    if(length(v) != 3) stop(nm, " priors need c(median, log.sd[, on])", call. = FALSE); v }))
  x <- matrix(as.numeric(x), ncol = 3, dimnames = list(rownames(x), c("median", "log.sd", "on")))
  if(any(x[, "on"] == 1 & !(x[, "median"] > 0 & x[, "log.sd"] > 0))) stop(nm, " priors need median > 0 and log.sd > 0", call. = FALSE)
  x
}

# maximum catch for the K prior and the K.max cap: raw maximum (JABBA) or maximum of the kernel-smoothed catch (CMSY, bw 3)
.max_catch <- function(data, type = "max"){
  C <- data$catch; ok <- is.finite(C)
  if(!any(ok)) return(NA_real_)
  if(type == "max" || sum(ok) < 4) return(max(C[ok]))
  y <- data$years[ok]
  max(stats::ksmooth(y, C[ok], kernel = "normal", bandwidth = 3, n.points = length(y))$y)
}

#' Predicted K / max catch from r or the generation time
#'
#' Median and log-SD of K / max catch predicted from the intrinsic rate of increase r or the generation time GT, from the
#' between-stock relations in 73 ICES stocks (K = unfished exploitable biomass, Beverton-Holt equilibrium):
#' K / max catch = 3.04 r^-0.583 (residual log-SD 0.55) and 1.71 GT^0.58 (0.57). Without r or GT, the empirical median
#' 5.9 (log-SD 0.55). The median is the centre of a lognormal; for an expected value use median x exp(log.sd^2/2).
#' Use with \code{jb3_priors(K = c(median x max catch, log.sd))} (two-sided) or \code{K.max = c(median, log.sd)} (one-sided,
#' as \code{K = "half+r"}).
#' @param r median of the r prior
#' @param GT generation time (years), used when \code{r} is NULL
#' @param maxC optional maximum catch: then also returns the median of K in biomass
#' @return named vector: median (multiple of max catch), log.sd, and K (if maxC given)
#' @examples
#' jb3_kr(r = 0.2)
#' jb3_kr(GT = 10, maxC = 1000)
#' @export
jb3_kr <- function(r = NULL, GT = NULL, maxC = NULL){
  v <- if(!is.null(r)) c(median = 3.04*r^-0.583, log.sd = 0.55) else if(!is.null(GT)) c(median = 1.71*GT^0.58, log.sd = 0.57) else c(median = 5.9, log.sd = 0.55)
  if(!is.null(maxC)) v <- c(v, K = v[["median"]]*maxC)
  v
}
