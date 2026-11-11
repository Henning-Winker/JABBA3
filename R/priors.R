#' Priors for JABBA3
#'
#' One named list of priors, each \code{c(median, log.sd, on)}: lognormal with the given median and log-scale sd, and an
#' on/off switch (\code{on = 0} gives a flat prior on the log scale). The defaults follow JABBA
#' (\code{build_jabba()}), except the process error prior.
#'
#' @param data optional \code{jb3data}, used for the default K prior (8 x max catch, CV 1, as JABBA)
#' @param r intrinsic growth rate, on by default
#' @param K carrying capacity; \code{NULL} = 8 x max catch with log.sd 0.833 (JABBA: CV 1)
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
#' @param sdf SD of the random walk in log harvest rate, used for years without catch (\code{jb3_data(catch.gaps =
#'   "estimate")}): lognormal prior c(median, log.sd, on), default median 0.3 (log.sd 0.5). SPiCT benchmarks often use a loose
#'   prior with the catch SD fixed, e.g. c(4, 0.5, 1), together with \code{jb3_model(frw.all = TRUE)}.
#' @param creep normal prior c(mean, sd) on the annual effort creep, used with \code{jb3_model(effort.creep = "estimate")}
#' @param psi alias of \code{bkinit} (JABBA name)
#' @param sdc catch SD prior c(median, log.sd, on) when \code{jb3_model(sdc.est = TRUE)} (off by default: flat on log sdc)
#' @param K.max optional soft upper prior on K, to avoid implausibly large biomass when the data say little about
#'   the scale: no penalty below the cap, and a half-normal on log K above it. \code{"r"} sets the cap from r (the median of
#'   the r prior) at the 95th percentile of K / max catch across 72 ICES stocks (K = unfished exploitable biomass):
#'   cap = 7.6 x max catch x r^-0.587 (about 19 x max catch at r = 0.2, 12 x at 0.4, 7.6 x at 1); a number gives the cap
#'   as a multiple of the maximum catch; \code{c(multiple, log.sd)} also sets the SD of the half-normal (default 0.3).
#'   \code{NULL} (default) = no cap
#' @return list of class \code{jb3priors}
#' @examples
#' jb3_priors(r = c(0.3, 0.3))
#' @export
jb3_priors <- function(data = NULL, r = c(0.2, 0.5, 1), K = NULL, bkinit = c(0.9, 0.246, 1), sdb = "fishlife", sigR = 0.5,
                       q = NULL, m = c(NA, 0.3, 0), kdyn = NULL, rho = c(0.7, 0.15, 0), sdi = c(0.2, 0.5, 0), status = NULL, psi = NULL,
                       creep = c(0.02, 0.02), sdf = c(0.3, 0.5, 1), K.max = NULL, sdc = c(0.1, 1, 0)){
  args <- .fn_args()
  pr3 <- function(x, nm){
    x <- as.numeric(x); if(length(x) == 2) x <- c(x, 1)
    if(length(x) != 3) stop("prior ", nm, " must be c(median, log.sd[, on])", call. = FALSE)
    stats::setNames(x, c("median", "log.sd", "on"))
  }
  if(!is.null(psi)) bkinit <- psi
  r <- pr3(r, "r")
  if(is.null(K)) K <- if(!is.null(data)) c(8*max(data$catch, na.rm = TRUE), cv2sd(1), 1) else c(NA, cv2sd(1), 1)
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
  K.max <- if(is.null(K.max)) c(mult = NA, log.sd = NA, on = 0) else if(identical(K.max, "r"))
    c(mult = 7.6*r[["median"]]^-0.587, log.sd = 0.3, on = 1) else {
      if(!is.numeric(K.max) || K.max[1] <= 0) stop("K.max must be \"r\", a multiple of the maximum catch, or c(multiple, log.sd)", call. = FALSE)
      c(mult = K.max[1], log.sd = if(length(K.max) > 1) K.max[2] else 0.3, on = 1) }
  structure(list(r = r, K = K, bkinit = bkinit, sdb = sdb, sdb.type = sdb.type, sigR = sigR, q = q, m = m, kdyn = kdyn, rho = rho,
                 sdi = sdi, sdf = sdf, sdc = stats::setNames(as.numeric(sdc), c("median", "log.sd", "on")), status = status, creep = stats::setNames(as.numeric(creep[1:2]), c("mean", "sd")), K.max = K.max),
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
