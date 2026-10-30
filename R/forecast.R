# Forecasts: forward control (fwdctrl), projections, weighted ensembles and the hockey-stick harvest control rule

# ---------------------------------------------------------------------------------------------------------------------
# core: stochastic forward projection of one fit's draws under a control data.frame (one row per projection year)
# Timing: the first transition (year N -> N+1) uses the estimated catch of the last assessment year N; the control applies
# from year N+1, with C_y = H_y * EB_y on the start-of-year biomass of y (the HCR uses B_y/Bmsy).
# ---------------------------------------------------------------------------------------------------------------------
.hockey <- function(b, ftgt = 1, btrigger = 0.5, bclose = 0, fmin = 0.001)
  ifelse(b <= bclose, fmin, ifelse(b < btrigger, (b - bclose)*(ftgt - fmin)/(btrigger - bclose) + fmin, ftgt))

.fwd_core <- function(fit, dr, ctrl, procerr = "ar1", hcr = list()){
  d <- fit$setup$dat; N <- d$N; nd <- nrow(dr$B); ny <- nrow(ctrl)
  pp <- dr$pars; K <- pp$K; r <- pp$r; m <- pp$m; BmsyK <- pp$Bmsy/pp$K; Hmsy <- pp$Fmsy; sdb <- pp$sdb; kd <- pp$kdyn
  rho <- fit$report$rho %||% 0; D <- dr$`Process deviation`[, N]
  rem <- function(p) if(d$eb.rem) .rho_curve(p, d$ebm[d$ebcc, ]) else 1
  sdf <- function(x) if(d$dyn) sdb*pmin(pmax(x, d$xlim[1]), d$xlim[2])^(-kd) else sdb
  hc <- utils::modifyList(list(btrigger = 0.5, bclose = 0, fmin = 0.001), hcr)
  P <- cbind(dr$`B/K`, matrix(NA, nd, ny)); Cm <- Hm <- matrix(NA, nd, ny)
  Cprev <- dr$Catch[, N]                                                       # last assessment year's catch
  for(t in seq_len(ny)){
    tt <- N + t; Pp <- P[, tt - 1]; Pl <- P[, max(1, tt - 1 - d$lag0)]
    pm <- log(pmax(Pp + r/(m - 1)*Pl*(1 - Pl^(m - 1)) - Cprev/(K*rem(Pp)), 0.001))
    v <- sdf(Pp/BmsyK)^2
    D <- switch(procerr,
      ar1 = if(isTRUE(d$rw)) stats::rnorm(nd, D, d$rw.scale*sqrt(v)) else stats::rnorm(nd, rho*D, sqrt(v*(1 - rho^2))),
      mean = if(isTRUE(d$rw)) D else rho*D,
      none = rep(0, nd))
    if(procerr != "ar1") v <- 0*v
    P[, tt] <- exp(pm - (if(d$biascor) 0.5*v else 0) + D)
    # control in year N + t, on the start-of-year biomass
    Pn <- P[, tt]; avail <- Pn*K*rem(Pn); q <- ctrl$quant[t]; val <- ctrl$value[t]
    H <- switch(q,
      Catch = pmin(val, 0.95*avail)/avail,
      F = pmin(rep(val, nd), 0.95),
      Fmsy = pmin(val*Hmsy, 0.95),
      HCR = pmin(.hockey(Pn/BmsyK, val, hc$btrigger, hc$bclose, hc$fmin)*Hmsy, 0.95),
      stop("unknown quant: ", q, call. = FALSE))
    Cm[, t] <- H*avail; Hm[, t] <- H/Hmsy; Cprev <- Cm[, t]
  }
  j <- N + seq_len(ny)
  list(`B/Bmsy` = P[, j, drop = FALSE]/BmsyK, B = P[, j, drop = FALSE]*K, Catch = Cm, `F/Fmsy` = Hm, F = Hm*Hmsy,
       years = max(fit$years) + seq_len(ny))
}

# members of a fit or ensemble: list(fit, draws)
.fwd_members <- function(x, ndraws){
  if(inherits(x, "jb3ensemble")) return(lapply(seq_along(x$fits), function(i) list(fit = x$fits[[i]], dr = x$memberdraws[[i]])))
  list(list(fit = x, dr = x$draws %||% jb3_draws(x, ndraws)))
}

.rbind_iters <- function(lst){
  out <- lapply(c("B/Bmsy", "B", "Catch", "F/Fmsy", "F"), function(q) do.call(rbind, lapply(lst, `[[`, q)))
  names(out) <- c("B/Bmsy", "B", "Catch", "F/Fmsy", "F"); out$years <- lst[[1]]$years; out
}

# status quo from the (pooled) draws: mean of the last nsq years
.sq <- function(x, nsq){
  dr <- if(inherits(x, "jb3ensemble")) x$draws else x$draws %||% NULL
  N <- length(x$years); j <- seq(N - nsq + 1, N)
  csq <- mean(utils::tail(x$data$catch, nsq))
  fsq <- if(!is.null(dr)) stats::median(rowMeans(dr$`F/Fmsy`[, j, drop = FALSE]))
         else mean(x$ts$est[x$ts$quantity == "F/Fmsy"][j])
  c(Csq = csq, Fsq = fsq)
}

.check_ctrl <- function(ctrl, x, nsq){
  ctrl <- as.data.frame(ctrl, stringsAsFactors = FALSE)
  if(!all(c("value", "quant") %in% names(ctrl))) stop("fwdctrl needs columns value and quant (and year)", call. = FALSE)
  y0 <- max(x$years)
  if(is.null(ctrl$year)) ctrl$year <- y0 + seq_len(nrow(ctrl))
  ctrl <- ctrl[order(ctrl$year), ]
  if(ctrl$year[1] != y0 + 1 || any(diff(ctrl$year) != 1))
    stop("fwdctrl years must be consecutive from ", y0 + 1, " (the year after the last assessment year)", call. = FALSE)
  sq <- .sq(x, nsq)
  i <- ctrl$quant == "Csq"; ctrl$value[i] <- ctrl$value[i]*sq["Csq"]; ctrl$quant[i] <- "Catch"
  i <- ctrl$quant == "Fsq"; ctrl$value[i] <- ctrl$value[i]*sq["Fsq"]; ctrl$quant[i] <- "Fmsy"
  i <- ctrl$quant %in% c("F/Fmsy", "Ftgt"); ctrl$quant[i] <- "Fmsy"
  ctrl
}

.ctrl_label <- function(ctrl){
  q <- utils::tail(ctrl$quant, 1); v <- utils::tail(ctrl$value, 1)
  switch(q, Catch = paste("Catch =", signif(v, 3)), F = paste("F =", signif(v, 3)), Fmsy = paste("F/Fmsy =", signif(v, 3)),
         HCR = paste("HCR Ftgt =", signif(v, 3)), paste(q, "=", signif(v, 3)))
}

#' Forward projections with a forecast control (fwdctrl)
#'
#' The JABBA3 equivalent of JABBA's \code{fw_jabba_ctrl()}: stochastic projections of a fit or a weighted ensemble under one
#' or more forecast controls. A control is a \code{data.frame(year, value, quant)} with one row per projection year,
#' starting in the year after the last assessment year; \code{quant} is one of
#' \describe{
#'   \item{"Catch"}{catch (tonnes)}
#'   \item{"Fmsy"}{F relative to Fmsy (alias "F/Fmsy")}
#'   \item{"F"}{harvest rate}
#'   \item{"Csq", "Fsq"}{multiplier of the status quo catch or F/Fmsy (mean of the last \code{nsq} years)}
#'   \item{"HCR"}{the hockey-stick rule (\code{\link{jb3_hcr}()}) applied draw by draw to the start-of-year B/Bmsy, with
#'     \code{value} = Ftgt/Fmsy and the rule settings in \code{hcr}}
#' }
#' The first transition (last assessment year to the next) uses the estimated catch of the last assessment year, so the
#' biomass in the first projection year is a one-step-ahead prediction. Process error follows the fitted model
#' (\code{procerr}, see \code{\link{jb3_project}()}).
#' @param x a \code{jabba3} fit or a \code{\link{jb3_ensemble}()}
#' @param fwdctrl a control data.frame or a (named) list of them, one per scenario
#' @param ndraws draws when the fit has none
#' @param procerr "ar1" (default), "mean" or "none"
#' @param nsq years for the status quo (Csq, Fsq)
#' @param hcr list of hockey-stick settings for quant = "HCR": btrigger, bclose, fmin (as fractions of Bmsy and Fmsy)
#' @param seed random seed
#' @return a \code{jb3proj} object (as \code{jb3_project()}), with the controls in \code{ctrl}
#' @examples
#' \dontrun{
#' csq <- mean(tail(bet$catch$catch, 3))
#' fc <- list(Fmsy = data.frame(year = 2018:2021, value = c(csq, csq, 1, 1), quant = c("Catch", "Catch", "Fmsy", "Fmsy")),
#'            Csq  = data.frame(year = 2018:2021, value = 1, quant = "Csq"))
#' stf <- jb3_fwd(fit, fc)
#' jbplot_projection(stf); jb3_fwdtable(stf)
#' }
#' @export
jb3_fwd <- function(x, fwdctrl, ndraws = 1000, procerr = c("ar1", "mean", "none"), nsq = 3, hcr = list(), seed = 123){
  procerr <- match.arg(procerr)
  if(is.data.frame(fwdctrl)) fwdctrl <- list(fwdctrl)
  ctrl <- lapply(fwdctrl, .check_ctrl, x = x, nsq = nsq)
  lab <- names(ctrl); if(is.null(lab)) lab <- rep("", length(ctrl))
  lab[lab == ""] <- vapply(ctrl[lab == ""], .ctrl_label, "")
  lab <- make.unique(lab); names(ctrl) <- lab
  mem <- .fwd_members(x, ndraws)
  set.seed(seed)
  out <- lapply(ctrl, function(cc) .rbind_iters(lapply(mem, function(mm) .fwd_core(mm$fit, mm$dr, cc, procerr, hcr))))
  last <- lapply(ctrl, function(cc) cc[nrow(cc), ])
  sc <- data.frame(type = vapply(last, function(z) switch(z$quant, Fmsy = "F/Fmsy", z$quant), ""),
                   value = vapply(last, function(z) z$value, 0), label = lab, row.names = NULL, stringsAsFactors = FALSE)
  # TAC year: first year in which a control departs from its first-year (intermediate) setting
  tac <- vapply(ctrl, function(cc){ k <- which(cc$quant != cc$quant[1] | cc$value != cc$value[1]); if(length(k)) cc$year[k[1]] else NA }, 0)
  tac <- if(all(is.na(tac))) ctrl[[1]]$year[1] else min(tac, na.rm = TRUE)
  .proj_object(x, out, sc, ctrl = ctrl, tac.year = tac, procerr = procerr, hcr = hcr)
}

.proj_object <- function(x, out, sc, ...){
  N <- length(x$years)
  ts <- do.call(rbind, lapply(seq_along(out), function(s) do.call(rbind, lapply(c("B/Bmsy", "B", "Catch", "F/Fmsy", "F"), function(q){
    qs <- apply(out[[s]][[q]], 2, stats::quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
    data.frame(scenario = sc$label[s], year = out[[s]]$years, quantity = q, est = qs[2, ], lci = qs[1, ], uci = qs[3, ],
               p.above = if(q == "B/Bmsy") colMeans(out[[s]][[q]] >= 1) else if(q == "F/Fmsy") colMeans(out[[s]][[q]] > 1 + 1e-8) else NA)
  }))))
  dr <- if(inherits(x, "jb3ensemble")) x$draws else x$draws %||% NULL
  st <- if(!is.null(dr)) list(`B/Bmsy` = dr$`B/Bmsy`[, N], `F/Fmsy` = dr$`F/Fmsy`[, N], B = dr$`B/Bmsy`[, max(1, N - 1)])
        else list(`B/Bmsy` = x$ts$est[x$ts$quantity == "B/Bmsy"][N], `F/Fmsy` = x$ts$est[x$ts$quantity == "F/Fmsy"][N])
  st$year <- max(x$years)
  names(out) <- sc$label
  structure(c(list(ts = ts, iters = out, scenarios = sc, run = x$run, fit.ts = x$ts, refpts = x$refpts, status = st), list(...)),
            class = "jb3proj")
}

#' Projections
#'
#' Stochastic projections under constant F/Fmsy (default) or constant catch scenarios, from posterior draws (NUTS) or
#' pseudo-posterior draws at the mode (generated if the fit has none), with process error and the production lag. Works
#' for a single fit and for a weighted \code{\link{jb3_ensemble}()}. Optional intermediate years (\code{imp.years}) with
#' a fixed catch (default the status quo catch) or F/Fmsy precede the scenarios, as the intermediate years of short-term
#' forecasts. A wrapper around \code{\link{jb3_fwd}()}, which takes general forecast controls.
#'
#' The first transition uses the estimated catch of the last assessment year, so the biomass in the first projection year
#' is a one-step-ahead prediction; catches are capped at 95\% of the exploitable biomass.
#' @param fit a \code{jabba3} fit or a \code{jb3ensemble}
#' @param fmsy vector of F/Fmsy multipliers (one scenario each); default \code{c(0, 0.5, 0.75, 1, 1.25)}
#' @param catch alternatively, vector of constant catches
#' @param nyears number of projection years (including the intermediate years)
#' @param ndraws draws when the fit has none
#' @param procerr process error in the projection years, consistent with the fitted model:
#'   \code{"ar1"} (default) continues each draw's AR(1) process from its last estimated deviation,
#'   d_t+1 = rho d_t + e, e ~ N(0, sdb^2 (1 - rho^2)), so persistence carries into the forecast and the marginal SD stays sdb
#'   (with rho = 0 this is iid sampling as in JABBA's fw_jabba());
#'   \code{"mean"} uses the conditional mean only, d_t+h = rho^h d_T (no new innovations; uncertainty from the parameters
#'   and states only);
#'   \code{"none"} sets the deviations to zero (deterministic production; the bias correction is dropped too)
#' @param imp.years number of intermediate years before the scenarios start (default 0)
#' @param imp.catch catch in the intermediate years (default: status quo catch, mean of the last \code{nsq} years)
#' @param imp.fmsy alternatively, F/Fmsy in the intermediate years
#' @param nsq years for the status quo
#' @return list of class \code{jb3proj}: \code{ts} (scenario, year, quantity, quantiles), \code{iters}, \code{scenarios},
#'   \code{status}, \code{ctrl}, \code{tac.year}
#' @export
jb3_project <- function(fit, fmsy = NULL, catch = NULL, nyears = 10, ndraws = 1000, procerr = c("ar1", "mean", "none"),
                        imp.years = 0, imp.catch = NULL, imp.fmsy = NULL, nsq = 3){
  procerr <- match.arg(procerr)
  if(is.null(catch) && is.null(fmsy)) fmsy <- c(0, 0.5, 0.75, 1, 1.25)
  type <- if(!is.null(fmsy)) "Fmsy" else "Catch"; vals <- if(!is.null(fmsy)) fmsy else catch
  if(imp.years >= nyears) stop("imp.years must be smaller than nyears", call. = FALSE)
  imp <- if(imp.years > 0){
    if(!is.null(imp.fmsy)) data.frame(value = imp.fmsy, quant = "Fmsy")
    else if(!is.null(imp.catch)) data.frame(value = imp.catch, quant = "Catch") else data.frame(value = 1, quant = "Csq")
  }
  ctrl <- lapply(vals, function(v){
    cc <- data.frame(year = max(fit$years) + seq_len(nyears), value = v, quant = type, stringsAsFactors = FALSE)
    if(imp.years > 0){ cc$value[seq_len(imp.years)] <- imp$value; cc$quant[seq_len(imp.years)] <- imp$quant }
    cc })
  names(ctrl) <- paste(if(type == "Fmsy") "F/Fmsy" else "Catch", "=", signif(vals, 3))
  jb3_fwd(fit, ctrl, ndraws = ndraws, procerr = procerr, nsq = nsq)
}

# ---------------------------------------------------------------------------------------------------------------------
# Ensemble
# ---------------------------------------------------------------------------------------------------------------------

#' Weighted ensemble of fits
#'
#' Combines several fits (e.g. alternative production functions, priors, lag and AR(1) settings) into one weighted
#' ensemble, as the GFCM/ICCAT JABBA ensembles: each fit contributes a number of draws proportional to its weight (posterior
#' draws for NUTS fits, MVN draws at the mode otherwise), and the joint draws give the ensemble trajectories, reference
#' points and status. The result behaves like a fit for \code{jbplot_advice()}, \code{jbplot_trj()}, \code{jbplot_kobe()},
#' \code{jb3_FLStockR()}, \code{jb3_partable()}-style summaries and projections (\code{jb3_project()}, \code{jb3_fwd()},
#' \code{jb3_hcr()}), which run each member's own model on its share of the draws.
#' @param fits named list of \code{jabba3} fits (same years)
#' @param weights weights (default equal), e.g. from \code{\link{jb3_weights}()}
#' @param ndraws total number of draws
#' @param seed random seed
#' @return object of class \code{c("jb3ensemble", "jabba3")}
#' @examples
#' \dontrun{
#' runs <- list(Schaefer = f1, Fox = f2, Pella = f3)
#' ens <- jb3_ensemble(runs, weights = jb3_weights(runs)$weight)
#' jbplot_ensemble(ens); jbplot_advice(ens); jbplot_kobe(ens)
#' }
#' @export
jb3_ensemble <- function(fits, weights = NULL, ndraws = 3000, seed = 123){
  if(is.null(names(fits))) names(fits) <- vapply(fits, function(f) f$run, "")
  names(fits) <- make.unique(names(fits))
  yrs <- fits[[1]]$years
  if(!all(vapply(fits, function(f) identical(f$years, yrs), TRUE))) stop("all fits need the same years", call. = FALSE)
  w <- weights %||% rep(1, length(fits)); if(length(w) != length(fits) || any(w < 0) || sum(w) <= 0) stop("invalid weights", call. = FALSE)
  w <- w/sum(w); n <- round(ndraws*w)
  keep <- n > 0
  fits <- fits[keep]; w <- w[keep]; n <- n[keep]
  set.seed(seed)
  md <- lapply(seq_along(fits), function(i){ f <- fits[[i]]
    dr <- if(!is.null(f$draws) && nrow(f$draws$B) >= 1) f$draws else jb3_draws(f, n[i], seed = seed + i)
    k <- sample.int(nrow(dr$B), n[i], replace = nrow(dr$B) < n[i])
    .subset_draws(dr, k) })
  q <- c("B", "B/Bmsy", "F", "F/Fmsy", "B/K", "Catch", "Process deviation")
  dr <- lapply(stats::setNames(q, q), function(z) do.call(rbind, lapply(md, `[[`, z)))
  dr$pars <- do.call(rbind, lapply(md, `[[`, "pars")); dr$member <- rep(names(fits), n); dr$type <- "ensemble"
  ts <- do.call(rbind, lapply(q, function(z){ qs <- apply(dr[[z]], 2, stats::quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
    data.frame(year = yrs, quantity = z, est = qs[2, ], lci = qs[1, ], uci = qs[3, ], se.log = NA, run = NA_character_) }))
  med <- function(v) c(stats::median(v), stats::quantile(v, c(0.025, 0.975), names = FALSE))
  rp <- data.frame(par = c("Bmsy", "Fmsy", "MSY", "K", "BmsyK"), est = NA, lci = NA, uci = NA, se.log = NA, row.names = c("Bmsy", "Fmsy", "MSY", "K", "BmsyK"))
  for(nm in c("Bmsy", "Fmsy", "MSY", "K")) rp[nm, c("est", "lci", "uci")] <- med(dr$pars[[nm]])
  rp["BmsyK", c("est", "lci", "uci")] <- med(dr$pars$Bmsy/dr$pars$K)
  pn <- c("K", "r", "m", "sdb", "bkinit")
  pars <- data.frame(par = pn, est = NA, lci = NA, uci = NA, se.log = NA, row.names = pn)
  for(nm in pn) pars[nm, c("est", "lci", "uci")] <- med(dr$pars[[nm]])
  f1 <- fits[[1]]
  structure(list(stock = f1$stock, years = yrs, data = f1$data, ts = ts, pars = pars, refpts = rp, draws = dr, run = "Ensemble",
                 fits = fits, weights = stats::setNames(w, names(fits)), n = stats::setNames(n, names(fits)), memberdraws = md,
                 report = list(rho = NA), engine = "ensemble"),
            class = c("jb3ensemble", "jabba3"))
}

.subset_draws <- function(dr, k){
  out <- dr
  for(nm in names(dr)){ x <- dr[[nm]]
    if(is.matrix(x) || is.data.frame(x)) out[[nm]] <- x[k, , drop = FALSE] }
  out
}

#' @export
print.jb3ensemble <- function(x, ...){
  cat("JABBA3 ensemble of", length(x$fits), "runs,", sum(x$n), "draws\n")
  print(data.frame(run = names(x$fits), weight = round(x$weights, 3), draws = x$n, row.names = NULL), row.names = FALSE)
  N <- length(x$years)
  cat("  status", max(x$years), ": B/Bmsy", sprintf("%.3f", stats::median(x$draws$`B/Bmsy`[, N])),
      " F/Fmsy", sprintf("%.3f", stats::median(x$draws$`F/Fmsy`[, N])), "\n")
  invisible(x)
}

#' Ensemble weights from model diagnostics
#'
#' Scores each fit on diagnostic criteria, as used to weight JABBA ensembles: convergence (a fit that did not converge gets
#' weight 0), runs tests of the residuals (all indices pass), the RMSE of the residuals (below \code{rmse.max}), and
#' optionally the hindcast prediction skill (MASE < 1 for all indices, from \code{\link{jb3_hindcast}()}) and the
#' retrospective pattern (|Mohn's rho| of B/Bmsy within \code{rho.lim}, from \code{\link{jb3_retro}()}). The weight is the
#' proportion of criteria passed (normalised to sum to 1); \code{weight = "equal"} gives equal weights to all converged runs.
#' @param fits named list of fits
#' @param hindcast optional list of \code{jb3_hindcast()} outputs (same order)
#' @param retro optional list of \code{jb3_retro()} outputs (same order)
#' @param rmse.max RMSE threshold (log scale)
#' @param rho.lim Mohn's rho bounds for B/Bmsy
#' @param weight "score" or "equal"
#' @return data.frame(run, converged, runs, rmse, mase, mohn, score, weight)
#' @export
jb3_weights <- function(fits, hindcast = NULL, retro = NULL, rmse.max = 0.3, rho.lim = c(-0.15, 0.2), weight = c("score", "equal")){
  weight <- match.arg(weight)
  nm <- names(fits) %||% vapply(fits, function(f) f$run, "")
  out <- do.call(rbind, lapply(seq_along(fits), function(i){ f <- fits[[i]]
    conv <- f$convergence$tier %in% c("ok", "converged", "acceptable")
    dg <- tryCatch(jb3_diags(f, if(any(is.finite(f$index$osa))) "osa" else "ordinary"), error = function(e) NULL)
    runs <- if(is.null(dg)) NA else all(dg$runs.p > 0.05, na.rm = TRUE)
    rmse <- if(is.null(dg)) NA else sqrt(mean(f$index$residual^2, na.rm = TRUE)) < rmse.max
    mase <- if(!is.null(hindcast)) all(hindcast[[i]]$mase$MASE < 1, na.rm = TRUE) else NA
    mohn <- if(!is.null(retro)){ r <- retro[[i]]$rho$rho[retro[[i]]$rho$quantity == "B/Bmsy"]; r >= rho.lim[1] && r <= rho.lim[2] } else NA
    crit <- c(runs, rmse, mase, mohn)
    data.frame(run = nm[i], converged = conv, runs = runs, rmse = rmse, mase = mase, mohn = mohn,
               score = if(conv) mean(crit, na.rm = TRUE) else 0) }))
  out$weight <- if(weight == "equal") as.numeric(out$converged) else out$score
  if(sum(out$weight) > 0) out$weight <- out$weight/sum(out$weight)
  out
}

# ---------------------------------------------------------------------------------------------------------------------
# Hockey-stick harvest control rule
# ---------------------------------------------------------------------------------------------------------------------

#' Hockey-stick harvest control rule
#'
#' The GFCM/ICES-type hockey-stick rule (JABBA's \code{jb_hockeystick()}): F/Fmsy = \code{ftgt} above \code{btrigger}
#' (B/Bmsy), decreasing linearly to \code{fmin} at \code{bclose} (0 by default, as JABBA; set \code{bclose = blim} to close
#' the fishery below Blim). The rule is evaluated on the biomass in the first implementation (TAC) year, after the
#' intermediate years of a forecast control: first an \code{ftgt} forecast is run (or \code{x} is such a forecast), then
#' \itemize{
#'   \item \code{uncertainty = FALSE}: F_adv from the median B/Bmsy;
#'   \item \code{uncertainty = TRUE} (precautionary, the default): F_adv = mean over all draws of the rule applied to
#'     each draw's B/Bmsy.
#' }
#' The returned F_adv (F/Fmsy) can be used in a new control (\code{quant = "Fmsy"}); \code{quant = "HCR"} in
#' \code{\link{jb3_fwd}()} applies the rule draw by draw in every projection year instead.
#' @param x a fit or \code{jb3ensemble} (an Ftgt forecast is run with \code{imp.years} intermediate years at
#'   \code{imp.catch}, default status quo catch), or a \code{jb3proj} from \code{jb3_fwd()}
#' @param ftgt target F/Fmsy
#' @param btrigger B/Bmsy below which F is reduced (GFCM: Bthr = 0.5 Bmsy)
#' @param blim limit B/Bmsy (shown in the plot; default 0.3)
#' @param bclose B/Bmsy at which F reaches \code{fmin}
#' @param fmin minimum F/Fmsy
#' @param year evaluation year (default: the TAC year of the forecast)
#' @param uncertainty precautionary F_adv across draws
#' @param imp.years,imp.catch intermediate years and their catch when \code{x} is a fit
#' @param scenario scenario of \code{x} to use when \code{x} is a forecast (default the first)
#' @param ... passed to \code{jb3_fwd()} (e.g. \code{ndraws}, \code{procerr})
#' @return object of class \code{jb3hcr}: \code{Fadv} (F/Fmsy), \code{Fadv.det}, \code{Fadv.pa}, B/Bmsy draws in the
#'   evaluation year, the rule settings and the forecast
#' @export
jb3_hcr <- function(x, ftgt = 1, btrigger = 0.5, blim = 0.3, bclose = 0, fmin = 0.001, year = NULL, uncertainty = TRUE,
                    imp.years = 2, imp.catch = NULL, scenario = 1, ...){
  if(!inherits(x, "jb3proj")){
    cc <- data.frame(year = max(x$years) + seq_len(imp.years + 2), value = ftgt, quant = "Fmsy", stringsAsFactors = FALSE)
    if(imp.years > 0){ cc$value[seq_len(imp.years)] <- imp.catch %||% 1; cc$quant[seq_len(imp.years)] <- if(is.null(imp.catch)) "Csq" else "Catch" }
    x <- jb3_fwd(x, list(Ftgt = cc), ...)
  }
  it <- x$iters[[scenario]]; year <- year %||% x$tac.year %||% it$years[1]
  j <- match(year, it$years); if(is.na(j)) stop("year not in the forecast", call. = FALSE)
  b <- it$`B/Bmsy`[, j]
  fdet <- .hockey(stats::median(b), ftgt, btrigger, bclose, fmin); fpa <- mean(.hockey(b, ftgt, btrigger, bclose, fmin))
  structure(list(Fadv = if(uncertainty) fpa else fdet, Fadv.det = fdet, Fadv.pa = fpa, b = b, year = year,
                 rule = c(ftgt = ftgt, btrigger = btrigger, blim = blim, bclose = bclose, fmin = fmin),
                 status = x$status, fwd = x), class = "jb3hcr")
}

#' @export
print.jb3hcr <- function(x, ...){
  cat(sprintf("Hockey-stick HCR (Ftgt %.2f Fmsy, Btrigger %.2f Bmsy) in %d: B/Bmsy median %.3f, P(B < Btrigger) %.2f\n",
              x$rule["ftgt"], x$rule["btrigger"], x$year, stats::median(x$b), mean(x$b < x$rule["btrigger"])))
  cat(sprintf("  F_adv/Fmsy: %.3f (median B), %.3f (precautionary, across draws) -> %.3f\n", x$Fadv.det, x$Fadv.pa, x$Fadv))
  invisible(x)
}

#' Harvest control rule plot
#'
#' The hockey-stick rule in F/Fmsy against B/Bmsy (as \code{FLRef::plotGFCM()}), with Blim, Btrigger and Bmsy, the
#' distribution of B/Bmsy in the evaluation year (the TAC year after the intermediate years) and the advised F.
#' @param x a \code{jb3hcr} (from \code{jb3_hcr()}); or NULL to draw the rule only
#' @param ftgt,btrigger,blim,bclose,fmin rule settings when \code{x} is NULL
#' @param bmax maximum B/Bmsy shown
#' @export
jbplot_hcr <- function(x = NULL, ftgt = 1, btrigger = 0.5, blim = 0.3, bclose = 0, fmin = 0.001, bmax = 2){
  if(!is.null(x)){ r <- x$rule; ftgt <- r[["ftgt"]]; btrigger <- r[["btrigger"]]; blim <- r[["blim"]]; bclose <- r[["bclose"]]; fmin <- r[["fmin"]]
    bmax <- max(bmax, stats::quantile(x$b, 0.99)*1.05) }
  bb <- seq(0, bmax, length.out = 400); rule <- data.frame(b = bb, f = .hockey(bb, ftgt, btrigger, bclose, fmin))
  fmax <- ftgt*1.4
  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0, xmax = blim, ymin = 0, ymax = Inf, fill = "#de2d26", alpha = 0.25) +
    ggplot2::annotate("rect", xmin = blim, xmax = btrigger, ymin = 0, ymax = Inf, fill = "#fee391", alpha = 0.45) +
    ggplot2::annotate("rect", xmin = btrigger, xmax = Inf, ymin = 0, ymax = Inf, fill = "#2ca25f", alpha = 0.2) +
    ggplot2::geom_vline(xintercept = c(blim, btrigger, 1), linetype = 2, colour = c("red", "darkorange", "darkgreen")) +
    ggplot2::annotate("text", x = c(blim, btrigger, 1), y = fmax*0.97, label = c("Blim", "Btrigger", "Bmsy"), hjust = -0.1, size = 3,
                      colour = c("red", "darkorange", "darkgreen")) +
    ggplot2::geom_line(data = rule, ggplot2::aes(.data$b, .data$f), linewidth = 1.1)
  if(!is.null(x)){
    q <- stats::quantile(x$b, c(0.05, 0.5, 0.95))
    dens <- stats::density(x$b, from = 0, to = bmax); dd <- data.frame(b = dens$x, y = dens$y/max(dens$y)*0.25*ftgt)
    p <- p + ggplot2::geom_area(data = dd, ggplot2::aes(.data$b, .data$y), fill = "dodgerblue", alpha = 0.35) +
      ggplot2::annotate("segment", x = q[1], xend = q[3], y = x$Fadv, yend = x$Fadv, colour = "dodgerblue4", linewidth = 0.8) +
      ggplot2::annotate("point", x = q[2], y = x$Fadv, shape = 21, fill = "white", colour = "dodgerblue4", size = 3.5, stroke = 1.2) +
      ggplot2::annotate("label", x = q[2], y = x$Fadv, label = sprintf("F_adv = %.2f Fmsy (%d)", x$Fadv, x$year), vjust = -0.8,
                        size = 3, label.size = 0, fill = "white")
  }
  p + ggplot2::coord_cartesian(xlim = c(0, bmax), ylim = c(0, fmax), expand = FALSE) +
    ggplot2::labs(x = expression(B/B[MSY]), y = expression(F/F[MSY]),
                  subtitle = if(!is.null(x)) "B/Bmsy distribution in the TAC year (blue), median and 90% interval" else NULL) +
    theme_jabba3()
}

# ---------------------------------------------------------------------------------------------------------------------
# Forecast table (GFCM fishing opportunities; as FLRef::fwd2type1)
# ---------------------------------------------------------------------------------------------------------------------

#' Fishing opportunities table
#'
#' The GFCM/ICES-style summary of fishing opportunities by forecast scenario (as \code{FLRef::fwd2type1()} used in the
#' GFCM JABBA ensemble guides): the catch and F in the evaluation (TAC) year, the biomass at the start of the following
#' year (\code{osa = TRUE}, one-step-ahead response to the TAC), its change relative to the year before the TAC year, and
#' the probabilities (\%) of the biomass being above Blim, Bpa (= Bthr) and Btgt (= Bmsy).
#' @param proj a \code{jb3proj} from \code{jb3_fwd()} or \code{jb3_project()}
#' @param eval.yr evaluation (TAC) year; default the TAC year of the forecast control (or the first projection year)
#' @param osa evaluate the biomass one step ahead (year after the TAC year)
#' @param rel relative units (F/Ftgt, B/Btgt) as in the GFCM tables; FALSE gives catch, F and biomass in absolute units
#' @param bthr,blim Bpa (Bthr) and Blim as fractions of Bmsy
#' @param digits decimals
#' @return data.frame with columns Basis, C y, F/Ftgt y, B/Btgt y+1, B change, P(B>Blim), P(B>Bpa), P(B>Btgt)
#' @export
jb3_fwdtable <- function(proj, eval.yr = NULL, osa = TRUE, rel = TRUE, bthr = 0.5, blim = 0.3, digits = 3){
  py <- proj$iters[[1]]$years
  eval.yr <- eval.yr %||% proj$tac.year %||% py[1]
  evalB <- if(osa) eval.yr + 1 else eval.yr
  if(!evalB %in% py) stop("projection too short: biomass in ", evalB, " needed (increase nyears)", call. = FALSE)
  j <- match(eval.yr, py); jb <- match(evalB, py); j0 <- match(eval.yr - 1, py)
  df <- do.call(rbind, lapply(seq_along(proj$iters), function(s){ it <- proj$iters[[s]]
    bb <- it$`B/Bmsy`[, jb]; bref <- if(!is.na(j0)) it[[if(rel) "B/Bmsy" else "B"]][, j0] else
      if(rel) proj$status$`B/Bmsy` else proj$fit.ts$est[proj$fit.ts$quantity == "B"][length(unique(proj$fit.ts$year))]
    Bv <- if(rel) bb else it$B[, jb]
    data.frame(Basis = proj$scenarios$label[s],
               C = round(stats::median(it$Catch[, j]), max(0, digits - 3)),
               F = round(stats::median(if(rel) it$`F/Fmsy`[, j] else it$F[, j]), digits),
               B = round(stats::median(Bv), if(rel) digits else 0),
               dB = round(100*(stats::median(Bv)/stats::median(bref) - 1), 1),
               PBlim = round(100*mean(bb > blim), 1), PBpa = round(100*mean(bb > bthr), 1), PBtgt = round(100*mean(bb > 1), 1))
  }))
  names(df) <- c("Basis", paste0("C ", eval.yr), if(rel) paste0("F/Ftgt ", eval.yr) else paste0("F ", eval.yr),
                 if(rel) paste0("B/Btgt ", evalB) else paste0("B ", evalB), "B change", "P(B>Blim)", "P(B>Bpa)", "P(B>Btgt)")
  rownames(df) <- NULL
  df
}
