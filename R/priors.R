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
#' @param sdi index observation SD prior (only with \code{obs.error = "spict"}), off by default
#' @param status optional status priors: a data.frame(year, type, median, log.sd) with type one of "BK", "BBmsy",
#'   "FFmsy" (lower case and JABBA's "bk", "bbmsy", "ffmsy" accepted; year NA = terminal year of the fit), the output of
#'   \code{\link{jb3_bprior}()}, or a JABBA \code{b.prior} vector c(mean, cv, year, type) (the cv used as log-scale SD, as
#'   in JABBA; year NA = terminal year)
#' @param creep normal prior c(mean, sd) on the annual effort creep, used with \code{jb3_model(effort.creep = "estimate")}
#' @param psi alias of \code{bkinit} (JABBA name)
#' @return list of class \code{jb3priors}
#' @examples
#' jb3_priors(r = c(0.3, 0.3))
#' @export
jb3_priors <- function(data = NULL, r = c(0.2, 0.5, 1), K = NULL, bkinit = c(0.9, 0.246, 1), sdb = "fishlife", sigR = 0.5,
                       q = NULL, m = c(NA, 0.3, 0), kdyn = NULL, rho = c(0.7, 0.15, 0), sdi = c(0.2, 0.5, 0), status = NULL, psi = NULL,
                       creep = c(0.02, 0.02)){
  pr3 <- function(x, nm){
    x <- as.numeric(x); if(length(x) == 2) x <- c(x, 1)
    if(length(x) != 3) stop("prior ", nm, " must be c(median, log.sd[, on])", call. = FALSE)
    stats::setNames(x, c("median", "log.sd", "on"))
  }
  if(!is.null(psi)) bkinit <- psi
  r <- pr3(r, "r")
  if(is.null(K)) K <- if(!is.null(data)) c(8*max(data$catch), cv2sd(1), 1) else c(NA, cv2sd(1), 1)
  K <- pr3(K, "K"); bkinit <- pr3(bkinit, "bkinit")
  sdb.type <- if(is.character(sdb)) match.arg(sdb, c("fishlife", "igamma")) else "lognormal"
  sdb <- if(sdb.type == "fishlife") pr3(c(0.80*sigR*r[["median"]]^0.45, 0.3, 1), "sdb") else
    if(sdb.type == "igamma") c(median = NA, log.sd = NA, on = 1) else pr3(sdb, "sdb")
  if(!is.null(q)){ q <- matrix(q, ncol = 3, dimnames = list(NULL, c("median", "log.sd", "on"))) }
  m <- pr3(m, "m")
  if(is.null(kdyn)) kdyn <- c(if(r[["median"]] > 1.5) 0 else max(0, 0.064 - 0.155*log(r[["median"]])), 0.15, 1)
  kdyn <- stats::setNames(as.numeric(c(kdyn, 1)[1:3]), c("mean", "sd", "on"))
  sdi <- pr3(sdi, "sdi")
  rho <- stats::setNames(as.numeric(c(rho, 1)[1:3]), c("mean", "sd", "on"))
  status <- .status_df(status)
  structure(list(r = r, K = K, bkinit = bkinit, sdb = sdb, sdb.type = sdb.type, sigR = sigR, q = q, m = m, kdyn = kdyn, rho = rho,
                 sdi = sdi, status = status, creep = stats::setNames(as.numeric(creep[1:2]), c("mean", "sd"))), class = "jb3priors")
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
  tab <- rbind(r = x$r, K = x$K, bkinit = x$bkinit, sdb = x$sdb, m = x$m, sdi = x$sdi)
  cat("JABBA3 priors (lognormal: median, log.sd; on = 1)\n")
  print(signif(tab, 3))
  cat("sdb:", switch(x$sdb.type, fishlife = sprintf("FishLife centre 0.80 x sigR (%.2f) x r^0.45", x$sigR),
                     igamma = "diffuse inverse-gamma(0.001, 0.001) on the variance (JABBA legacy)", lognormal = "user lognormal"), "\n")
  cat("q:", if(is.null(x$q)) "flat on log q" else paste(nrow(x$q), "lognormal priors"), "\n")
  cat(sprintf("kdyn (dynamic sdb only): normal(%.2f, %.2f)%s\n", x$kdyn[1], x$kdyn[2], if(x$kdyn[3] == 0) " fixed at mean" else ""))
  if(!is.null(x$status)){ cat("status priors:\n"); print(x$status) }
  invisible(x)
}
