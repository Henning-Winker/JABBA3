# Convert a JABBA (JAGS) fit into a jabba3 fit object without refitting

#' Convert a JABBA fit into a JABBA3 fit object (no refit)
#'
#' Populates a \code{jabba3} fit object with the results of a JABBA fit (\code{JABBA::fit_jabba()}), so that the JABBA3
#' plots, tables and summaries can be used on it and the fit can be refitted with JABBA3 by \code{update()}. Data,
#' priors and model settings are reconstructed with \code{\link{jbinput2jabba3}()}. The estimates are the JABBA
#' posterior: trajectories (B, F, B/Bmsy, F/Fmsy, B/K, catch, process deviations) as posterior draws (from the JABBA
#' \code{kbtrj} element, as \code{fit$draws}, \code{type = "jags"}) and their medians and 95\% intervals (\code{fit$ts}),
#' parameters and reference points as posterior medians and intervals, and the index fits, residuals and posterior
#' predictive intervals as reported by JABBA. Nothing is estimated by JABBA3, and the fit has no RTMB object.
#'
#' \code{update(x)} refits the same data, priors and model with JABBA3 (\code{engine = "rtmb"}, the posterior mode; give
#' \code{engine = "nuts"} for the posterior), with the process-error prior proxy of \code{jbinput2jabba3()} when JABBA
#' used the diffuse inverse-gamma prior; any \code{jabba3()} argument and any prior or model change can be passed to
#' \code{update()} as usual. Functions that need the RTMB object (projections, retrospectives, hindcasts, draws, OSA
#' residuals) work on the refitted object.
#' @param jbfit output of \code{JABBA::fit_jabba()}
#' @param run run name (default: the JABBA scenario)
#' @param ... passed to \code{\link{jbinput2jabba3}()} (e.g. \code{igamma.proxy}, \code{sdb.inflate}, \code{q.uniform})
#' @return a \code{jabba3} fit with \code{engine = "jags"} and the original JABBA results in \code{$jabba}
#' @examples
#' \dontrun{
#' jb   <- JABBA::build_jabba(catch = bet$catch, cpue = bet$cpue, se = bet$se, model.type = "Fox")
#' jbf  <- JABBA::fit_jabba(jb)
#' f.jb <- jabba2jabba3(jbf)              # JABBA results in a jabba3 object: JABBA3 plots
#' jbplot_trj(f.jb); jbplot_kobe(f.jb); jbplot_fits(f.jb)
#' f.j3 <- update(f.jb)                   # refit with JABBA3
#' jbplot_trj(list(JABBA = f.jb, JABBA3 = f.j3))
#' }
#' @export
jabba2jabba3 <- function(jbfit, run = NULL, ...){
  if(is.null(jbfit$inputseries) || is.null(jbfit$settings)) stop("jbfit must be a JABBA::fit_jabba() object", call. = FALSE)
  x <- jbinput2jabba3(jbfit, ...)
  s <- jbfit$settings; yrs <- as.integer(x$data$years); N <- length(yrs)
  qfun <- function(v) c(stats::median(v, na.rm = TRUE), stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  sel <- function(v){ v <- v[is.finite(v) & v > 0]; if(length(v) > 1) stats::sd(log(v)) else NA_real_ }
  # posterior draws of the trajectories (kbtrj: one row per year and iteration)
  dr <- NULL; kb <- jbfit$kbtrj
  if(!is.null(kb)){
    if(!is.null(kb$type)) kb <- kb[kb$type == "fit", ]
    kb <- kb[kb$year %in% yrs, ]
    mat <- function(v){ w <- stats::reshape(kb[, c("iter", "year", v)], idvar = "iter", timevar = "year", direction = "wide")
      m <- as.matrix(w[, -1, drop = FALSE]); colnames(m) <- sub(paste0("^", v, "\\."), "", colnames(m)); m[, as.character(yrs), drop = FALSE] }
    pp <- jbfit$pars_posterior; rp <- jbfit$refpts_posterior; n <- length(unique(kb$iter))
    take <- function(v) if(!is.null(v) && length(v) >= n) v[seq_len(n)] else rep(NA_real_, n)
    qn <- grep("^q", names(pp), value = TRUE); tn <- grep("^tau2", names(pp), value = TRUE)
    dr <- list(B = mat("B"), `B/Bmsy` = mat("stock"), F = mat("H"), `F/Fmsy` = mat("harvest"), `B/K` = mat("BB0"),
               Catch = mat("Catch"), `Process deviation` = mat("Bdev"),
               pars = data.frame(K = take(pp$K), r = take(pp$r), bkinit = take(pp$psi), sdb = sqrt(take(pp$sigma2)),
                                 m = take(pp$m %||% rep(s$m %||% NA, n)), kdyn = NA_real_,
                                 Bmsy = take(rp$Bmsy), Fmsy = take(rp$Fmsy), MSY = take(rp$MSY)),
               q = if(length(qn)) as.matrix(pp[seq_len(n), qn, drop = FALSE]) else NULL,
               sdi = if(length(tn)) sqrt(as.matrix(pp[seq_len(n), tn, drop = FALSE])) else NULL, type = "jags")
  }
  # trajectories: posterior medians and 95% intervals
  tsq <- c(B = "B", F = "F", BBmsy = "B/Bmsy", FFmsy = "F/Fmsy", BB0 = "B/K")
  tr <- jbfit$timeseries; ty <- as.integer(dimnames(tr)[[1]]); it <- match(yrs, ty)
  mk <- function(q, est, lci, uci, se = NA_real_) data.frame(year = yrs, quant = q, est = est, lci = lci, uci = uci, se.log = se, run = NA_character_)
  ts <- if(!is.null(dr)) do.call(rbind, lapply(c("B", "F", "B/Bmsy", "F/Fmsy", "B/K", "Catch", "Process deviation"), function(q){
      v <- dr[[q]]; qs <- apply(v, 2, qfun); mk(q, qs[1, ], qs[2, ], qs[3, ], if(q == "Process deviation") apply(v, 2, stats::sd) else apply(v, 2, sel)) }))
    else rbind(do.call(rbind, lapply(names(tsq), function(k) mk(tsq[[k]], tr[it, "mu", k], tr[it, "lci", k], tr[it, "uci", k]))),
               mk("Catch", jbfit$est.catch$mu %||% x$data$catch, jbfit$est.catch$lci %||% NA, jbfit$est.catch$uci %||% NA),
               mk("Process deviation", log(tr[it, "mu", "procB"]), log(tr[it, "lci", "procB"]), log(tr[it, "uci", "procB"])))
  sdb.q <- qfun(sqrt(jbfit$pars_posterior$sigma2 %||% jbfit$pars["sigma2", "Median"]))
  sp <- if("SPt" %in% dimnames(tr)[[3]]) mk("Surplus production", tr[it, "mu", "SPt"], tr[it, "lci", "SPt"], tr[it, "uci", "SPt"]) else NULL
  ts <- rbind(ts, mk("Process SD", rep(sdb.q[1], N), rep(sdb.q[2], N), rep(sdb.q[3], N)), sp)
  rownames(ts) <- NULL
  # parameters and reference points (posterior medians, 95% intervals)
  pr <- function(nm, v) data.frame(par = nm, est = qfun(v)[1], lci = qfun(v)[2], uci = qfun(v)[3], se.log = sel(v), row.names = nm)
  pp <- jbfit$pars_posterior; P <- jbfit$pars
  pv <- function(nm, alt = NULL){ if(!is.null(pp[[nm]])) pp[[nm]] else if(nm %in% rownames(P)) c(P[nm, "LCI"], P[nm, "Median"], P[nm, "UCI"]) else alt }
  qn <- grep("^q", c(names(pp), rownames(P)), value = TRUE); qn <- unique(qn[qn != "qA"])
  tn <- unique(grep("^tau2", c(names(pp), rownames(P)), value = TRUE))
  pars <- rbind(pr("K", pv("K")), pr("r", pv("r")), pr("bkinit", pv("psi")), pr("sdb", sqrt(pv("sigma2"))),
                do.call(rbind, lapply(seq_along(qn), function(i) pr(paste0("q", i), pv(qn[i])))),
                if(length(tn)) do.call(rbind, lapply(seq_along(tn), function(i) pr(paste0("sdi", i), sqrt(pv(tn[i]))))),
                data.frame(par = "m", est = stats::median(pv("m", s$m %||% NA)), lci = NA, uci = NA, se.log = NA, row.names = "m"))
  rp <- jbfit$refpts_posterior
  refpts <- rbind(pr("Bmsy", rp$Bmsy %||% pv("SBmsy")), pr("Fmsy", rp$Fmsy %||% pv("Hmsy")), pr("MSY", rp$MSY %||% pv("MSY")),
                  pr("K", pv("K")), data.frame(par = "BmsyK", est = stats::median(rp$BmsyK %||% x$model$BmsyK), lci = NA, uci = NA, se.log = NA, row.names = "BmsyK"))
  # index fits (JABBA: posterior median and 95% interval of the fit, posterior predictive interval, total obs SD)
  ch <- jbfit$cpue.hat; pd <- jbfit$cpue.ppd; I <- as.matrix(x$data$I); nI <- ncol(I); inames <- colnames(I)
  hy <- as.integer(dimnames(ch)[[1]]); ih <- match(yrs, hy); sr <- jbfit$std.residuals
  index <- do.call(rbind, lapply(seq_len(nI), function(j){
    o <- I[, j]; h <- ch[ih, "mu", j]
    data.frame(index = inames[j], year = yrs, obs = o, hat = h, lci = ch[ih, "lci", j], uci = ch[ih, "uci", j], sd = ch[ih, "obserror", j],
               ppd_lci = if(!is.null(pd)) pd[ih, "lci", j] else NA, ppd_uci = if(!is.null(pd)) pd[ih, "uci", j] else NA,
               residual = log(o) - log(h), std.residual = if(!is.null(sr)) as.numeric(sr[j, match(yrs, as.integer(colnames(sr)))]) else NA, osa = NA_real_) }))
  index <- index[is.finite(index$obs), ]; rownames(index) <- NULL
  gw <- P[, "Geweke.p", drop = TRUE]; names(gw) <- rownames(P)
  conv <- list(code = 0, message = "JABBA (JAGS) fit: Geweke and Heidelberger-Welch diagnostics in $jabba$pars", maxgrad = NA, pdHess = NA,
               se.ok = TRUE, decrement = NA, tier = "jags", nll = NA, nstart = NA,
               geweke.min = suppressWarnings(min(gw, na.rm = TRUE)))
  setup <- tryCatch(jb3_setup(x$data, x$priors, x$model), error = function(e) NULL)
  fit <- list(stock = x$data$stock, years = yrs, data = x$data, ts = ts, pars = pars, refpts = refpts, convergence = conv,
              report = NULL, obj = NULL, opt = NULL, sdr = NULL, setup = setup, mucatch = setup$dat$mucatch %||% mean(x$data$catch, na.rm = TRUE),
              ts.season = NULL, nseason = 1L, run = run %||% jbfit$scenario %||% x$data$stock, priors = x$priors, model = x$model,
              engine = "jags", peels = 0, index.peels = 0, prior.only = FALSE, index = index, auxiliary = NULL, time = NA_real_,
              settings = list(engine = "rtmb", ndraws = 0, nuts = list(), osa = TRUE, workflow = "none", amat = NULL, se = TRUE,
                              peels = 0, index.peels = 0, nstart = 4, silent = TRUE, prior.only = FALSE),
              draws = dr, jabba = c(x$jabba, list(settings = s, stats = jbfit$stats, kobe = jbfit$kobe)))
  class(fit) <- c("jabba3", "list")
  fit
}
