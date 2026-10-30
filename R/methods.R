# Print, summary and FLR conversion

#' @export
print.jabba3 <- function(x, ...){
  cv <- x$convergence
  st <- x$ts[x$ts$quantity %in% c("B/Bmsy", "F/Fmsy") & x$ts$year == max(x$years), ]
  cat("JABBA3 fit:", x$run, "(", x$stock, ")  engine:", x$engine,
      if(!is.null(x$draws)) paste0(" + ", nrow(x$draws$B), " ", x$draws$type, " draws") else "", "\n")
  cat("  model:", x$model$shape, if(x$model$lag > 0) paste0("+ lag ", x$model$lag) else "",
      if(x$model$sdb.form == "dyn") "+ dynamic sdb" else "", if(!is.null(x$model$ebsb)) "+ Selex" else "", "\n")
  cat(sprintf("  convergence: %s (max gradient %.1e, Newton decrement %.1e)%s\n", cv$tier, cv$maxgrad, cv$decrement,
              if(!is.null(x$nuts)) sprintf("; NUTS Rhat max %.3f, ESS min %.0f, divergent %d", x$nuts$rhat.max, x$nuts$ess.min, x$nuts$divergent) else ""))
  p <- x$pars[, c("est", "lci", "uci")]; rp <- x$refpts[c("Bmsy", "Fmsy", "MSY"), c("est", "lci", "uci")]
  cat("  estimates (", if(x$engine == "nuts") "posterior median, 95% CrI" else "mode, 95% CI", "):\n", sep = "")
  print(signif(rbind(p, rp), 3))
  cat(sprintf("  status %d: B/Bmsy %.2f (%.2f-%.2f), F/Fmsy %.2f (%.2f-%.2f)\n", max(x$years),
              st$est[st$quantity == "B/Bmsy"], st$lci[st$quantity == "B/Bmsy"], st$uci[st$quantity == "B/Bmsy"],
              st$est[st$quantity == "F/Fmsy"], st$lci[st$quantity == "F/Fmsy"], st$uci[st$quantity == "F/Fmsy"]))
  if(!is.null(x$steps)) cat("  workflow steps:", paste(names(x$steps), collapse = " -> "), "\n")
  invisible(x)
}

#' @export
summary.jabba3 <- function(object, ...){
  out <- list(run = object$run, pars = object$pars, refpts = object$refpts, convergence = jb3_convergence(object),
              diags = if(!is.null(object$osa)) jb3_diags(object) else NULL,
              status = object$ts[object$ts$year == max(object$years) & object$ts$quantity %in% c("B", "F", "B/Bmsy", "F/Fmsy", "B/K"), ])
  class(out) <- "summary.jabba3"; out
}

#' @export
print.summary.jabba3 <- function(x, ...){
  cat("JABBA3 summary:", x$run, "\n\nConvergence (", attr(x$convergence, "tier"), ")\n", sep = ""); print(x$convergence[, 1:3], row.names = FALSE)
  cat("\nParameters\n"); print(signif(x$pars[, c("est", "lci", "uci")], 3))
  cat("\nReference points\n"); print(signif(x$refpts[, c("est", "lci", "uci")], 3))
  cat("\nTerminal status\n"); print(x$status[, c("quantity", "est", "lci", "uci")], row.names = FALSE, digits = 3)
  if(!is.null(x$diags)){ cat("\nResidual diagnostics (OSA)\n"); print(x$diags, row.names = FALSE, digits = 3) }
  invisible(x)
}

#' @export
print.jb3retro <- function(x, ...){ cat("JABBA3 retrospective,", length(x$fits) - 1, "peels. Mohn's rho:\n"); print(x$rho, row.names = FALSE, digits = 3); invisible(x) }
#' @export
print.jb3hindcast <- function(x, ...){ cat("JABBA3 hindcast,", x$n, "peels. MASE:\n"); print(x$mase, row.names = FALSE, digits = 3); invisible(x) }
#' @export
print.jb3proj <- function(x, ...){
  cat("JABBA3 projections:", nrow(x$scenarios), "scenarios\n")
  y <- max(x$ts$year); print(x$ts[x$ts$year == y & x$ts$quantity == "B/Bmsy", c("scenario", "year", "est", "lci", "uci", "p.above")], row.names = FALSE, digits = 3)
  invisible(x)
}

#' FLR outputs
#'
#' Trajectories as \code{FLQuants} (iterations = draws when the fit has them, otherwise the mode) and reference points as
#' \code{FLPar}, for use with ggplotFL, FLRef and mse.
#' @param fit a \code{jabba3} fit
#' @param draws use the draws (TRUE when available)
#' @return \code{jb3_flqs()}: FLQuants stock (B/Bmsy), harvest (F/Fmsy), biomass, F, catch; \code{jb3_flpar()}: FLPar
#' @export
jb3_flqs <- function(fit, draws = !is.null(fit$draws)){
  .need("FLCore", "FLR outputs need FLCore")
  yrs <- fit$years
  mk <- function(q){
    x <- if(draws) t(fit$draws[[q]]) else matrix(fit$ts$est[fit$ts$quantity == q], ncol = 1)
    FLCore::FLQuant(array(x, dim = c(1, length(yrs), 1, 1, 1, ncol(x))), dimnames = list(age = "all", year = yrs, iter = seq_len(ncol(x))))
  }
  FLCore::FLQuants(stock = mk("B/Bmsy"), harvest = mk("F/Fmsy"), biomass = mk("B"), F = mk("F"), catch = mk("Catch"))
}

#' @rdname jb3_flqs
#' @export
jb3_flpar <- function(fit){
  .need("FLCore", "FLR outputs need FLCore")
  rp <- fit$refpts; p <- fit$pars
  FLCore::FLPar(c(K = rp["K", "est"], r = p["r", "est"], m = p["m", "est"], Bmsy = rp["Bmsy", "est"], Fmsy = rp["Fmsy", "est"],
                  MSY = rp["MSY", "est"], BmsyK = rp["BmsyK", "est"], sdb = p["sdb", "est"],
                  B.Bmsy = utils::tail(fit$ts$est[fit$ts$quantity == "B/Bmsy"], 1),
                  F.Fmsy = utils::tail(fit$ts$est[fit$ts$quantity == "F/Fmsy"], 1)))
}

#' Convert a JABBA3 fit to an FLStockR
#'
#' Biomass-dynamic FLStockR (one "age"), as \code{FLRef::jabba2FLStockR()}: \code{stock} and \code{ssb()} are the biomass,
#' \code{catch} the (estimated) catch and \code{harvest} the harvest rate F = C/B; iterations are the posterior (NUTS) or
#' MVN draws, or the point estimates. The \code{refpts} slot holds Fmsy, Bmsy, MSY, Bthr, Blim, B0 (= K) and, with
#' \code{rel = TRUE}, everything is relative to Bmsy and Fmsy (Fmsy = Bmsy = 1). Use with \code{FLRef::plotAdvice()},
#' \code{ggplotFL} and \code{mse}.
#' @param fit a \code{jabba3} fit
#' @param bthr,blim biomass threshold and limit as fractions of Bmsy (defaults 0.5 and 0.3, SPiCT benchmark defaults)
#' @param rel relative units (B/Bmsy, F/Fmsy)
#' @param draws use the fit's draws as iterations (default if present); \code{ndraws} > 0 generates MVN draws otherwise
#' @param ndraws number of MVN draws to generate when the fit has none (0: point estimates)
#' @return FLStockR
#' @examples
#' \dontrun{
#' stk <- jb3_FLStockR(fit, ndraws = 500)
#' FLRef::plotAdvice(stk)
#' }
#' @export
jb3_FLStockR <- function(fit, bthr = 0.5, blim = 0.3, rel = FALSE, draws = !is.null(fit$draws), ndraws = 0){
  .need("FLCore", "jb3_FLStockR() needs FLCore")
  dr <- if(draws) fit$draws else if(ndraws > 0) jb3_draws(fit, ndraws) else NULL
  yrs <- fit$years; ny <- length(yrs)
  get <- function(q) if(!is.null(dr)) t(dr[[q]]) else matrix(fit$ts$est[fit$ts$quantity == q], ncol = 1)
  B <- get(if(rel) "B/Bmsy" else "B"); H <- get(if(rel) "F/Fmsy" else "F"); C <- get("Catch")
  ni <- ncol(B)
  flq <- function(x, v = NULL) FLCore::FLQuant(if(is.null(v)) array(x, dim = c(1, ny, 1, 1, 1, ni)) else v,
                                              dimnames = list(age = "1", year = yrs, iter = seq_len(ni)))
  one <- function(v) FLCore::FLQuant(v, dimnames = list(age = "1", year = yrs))
  stk <- FLCore::FLStock(stock.n = flq(B), catch.n = flq(C), landings.n = flq(C), discards.n = one(0),
                         stock.wt = one(1), catch.wt = one(1), landings.wt = one(1), discards.wt = one(1),
                         mat = one(1), m = one(1e-4), harvest = flq(H), m.spwn = one(0), harvest.spwn = one(0), name = fit$stock %||% "stock")
  FLCore::units(FLCore::harvest(stk)) <- "hr"
  FLCore::catch(stk) <- FLCore::computeCatch(stk); FLCore::landings(stk) <- FLCore::computeLandings(stk)
  FLCore::discards(stk) <- FLCore::computeDiscards(stk); FLCore::stock(stk) <- FLCore::computeStock(stk)
  stk <- FLCore::FLStockR(stk)
  rp <- fit$refpts
  bm <- if(rel) 1 else rp["Bmsy", "est"]; fm <- if(rel) 1 else rp["Fmsy", "est"]
  stk@refpts <- FLCore::FLPar(Fmsy = fm, Bmsy = bm, MSY = rp["MSY", "est"], Bthr = bthr*bm, Blim = blim*bm,
                              B0 = if(rel) rp["K", "est"]/rp["Bmsy", "est"] else rp["K", "est"])
  stk
}
