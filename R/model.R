#' Model structure for JABBA3
#'
#' Structural settings, separate from the priors (as \code{settings} in FLicc).
#'
#' @param shape "pella" (shape from \code{BmsyK}), "schaefer" (Bmsy/K = 0.5) or "fox" (m = 1.001)
#' @param BmsyK Bmsy/K for the Pella-Tomlinson form (0.4 default, as JABBA)
#' @param lag production lag in years (Tier 1: age at maturity); production in year t uses the depletion of t - lag
#' @param sdb.form "const" (default) or "dyn": process error SD sdb_t = sdb x (B/Bmsy)^-k (Tier 2)
#' @param sdb.fixed fixed process error SD at Bmsy instead of estimating it (Tier 2 sets the tuned level here)
#' @param xlim bounds of B/Bmsy in the dynamic SD
#' @param rho AR(1) coefficient of the process deviations (0 = independent, default). sdb stays the marginal SD:
#'   d_t ~ N(rho d_{t-1}, sdb_t^2 (1 - rho^2)). Estimated instead when \code{priors$rho} is switched on.
#' @param rw random walk in the process deviations (the conditional rho = 1 limit): d_t ~ N(d_{t-1}, (rw.scale sdb_t)^2)
#' @param rw.scale innovation SD of the random walk relative to sdb (e.g. sqrt(1 - 0.7^2) to match AR(1) rho = 0.7)
#' @param biascor lognormal process error bias correction (log mean - var/2)
#' @param start "all" (default): process error in all years; "index": deterministic before the first index year
#'   (JABBA \code{proc.dev.all = FALSE})
#' @param obs.error "additive" (default, JABBA): total SD = sqrt(se^2 + sdi^2); "spict": SD = sdi x relative se weights
#' @param se.min fixed minimum observation SD; NULL = JABBA default (0.1 without SEs, 0.001 with SEs)
#' @param se.na SE used where an index has SEs but a year has none (JABBA: 0.3)
#' @param sets.q,sets.var catchability and observation-variance groups (one entry per index; default one each)
#' @param catch.error estimate catch with lognormal error (random effects with SD \code{sdc}); on by default, as JABBA
#'   (\code{add.catch.CV = TRUE})
#' @param sdc catch observation SD (log scale) when \code{catch.error = TRUE}: a constant or a vector along the years
#'   (year-specific catch uncertainty; see also \code{jb3_data(catch.se =)}, which takes precedence where given)
#' @param ebsb Selex: EB/SB ratio curve c(u1, u2, u3, P1, P2), a named list of curves (one per selectivity) or the output
#'   of \code{\link{jb3_selex}()}; the state
#'   is then spawning biomass, indices are modelled as q x SB x (EB/SB)(SB/SB0)
#' @param ebsb.index curve per index (number or name; 0 = spawning-biomass index); default curve 1 for all
#' @param ebsb.catch curve for the removals (default 1)
#' @param ebsb.removal catch removes C / (EB/SB) from spawning biomass (default TRUE)
#' @param q.uniform legacy JABBA uniform prior on q (density proportional to q on log q)
#' @param P.bound,K.bound soft bounds on B/K and K (JABBA defaults)
#' @param sdi.bound soft upper bound on sdi (JABBA: 1)
#' @return list of class \code{jb3model}
#' @export
jb3_model <- function(shape = c("pella", "schaefer", "fox"), BmsyK = 0.4, lag = 0, sdb.form = c("const", "dyn"),
                      sdb.fixed = NULL, xlim = c(0.2, 3), rho = 0, rw = FALSE, rw.scale = 1, biascor = FALSE, start = c("all", "index"),
                      obs.error = c("additive", "spict"), se.min = NULL, se.na = 0.3, sets.q = NULL, sets.var = NULL,
                      catch.error = TRUE, sdc = 0.1, ebsb = NULL, ebsb.index = NULL, ebsb.catch = 1, ebsb.removal = TRUE,
                      q.uniform = FALSE, P.bound = c(0.02, 1.3), K.bound = c(0.01, 1e10), sdi.bound = 1){
  shape <- match.arg(shape); sdb.form <- match.arg(sdb.form); start <- match.arg(start); obs.error <- match.arg(obs.error)
  if(shape == "schaefer") BmsyK <- 0.5
  m <- switch(shape, schaefer = 2, fox = 1.001, pella = .m_from_bmsyk(BmsyK))
  if(shape == "fox") BmsyK <- m^(-1/(m - 1))
  lag <- as.integer(round(lag)); if(lag < 0) stop("lag must be >= 0", call. = FALSE)
  if(inherits(ebsb, "jb3selex")) ebsb <- ebsb$ebsb                       # from jb3_selex()
  if(!is.null(ebsb)){
    if(!is.list(ebsb)) ebsb <- list(ebsb)
    if(is.data.frame(ebsb)) ebsb <- lapply(seq_len(nrow(ebsb)), function(i) unlist(ebsb[i, ]))
    if(is.null(names(ebsb))) names(ebsb) <- paste0("sel", seq_along(ebsb))
    ebsb <- lapply(ebsb, function(e){ e <- as.numeric(e); if(length(e) != 5) stop("each ebsb curve is c(u1, u2, u3, P1, P2)", call. = FALSE); e })
  }
  if(rho < 0 || rho >= 1) stop("rho must be in [0, 1)", call. = FALSE)
  structure(list(shape = shape, rho = rho, rw = rw, rw.scale = rw.scale, BmsyK = BmsyK, m = m, lag = lag, sdb.form = sdb.form, sdb.fixed = sdb.fixed, xlim = xlim,
                 biascor = biascor, start = start, obs.error = obs.error, se.min = se.min, se.na = se.na, sets.q = sets.q,
                 sets.var = sets.var, catch.error = catch.error, sdc = sdc, ebsb = ebsb, ebsb.index = ebsb.index,
                 ebsb.catch = ebsb.catch, ebsb.removal = ebsb.removal, q.uniform = q.uniform, P.bound = P.bound,
                 K.bound = K.bound, sdi.bound = sdi.bound), class = "jb3model")
}

#' @export
print.jb3model <- function(x, ...){
  cat("JABBA3 model:", x$shape, sprintf("(m = %.3f, Bmsy/K = %.3f)", x$m, x$BmsyK),
      "\n  production lag:", x$lag, "\n  process error:", if(x$sdb.form == "const") "constant sdb" else "dynamic sdb x (B/Bmsy)^-k",
      if(!is.null(x$sdb.fixed)) sprintf("(level fixed at %.3f)", x$sdb.fixed) else "", if(x$biascor) "+ bias correction" else "", if(x$rho > 0) sprintf("+ AR(1) rho = %.2f", x$rho) else "",
      if(x$start == "index") "\n  deterministic before the first index year" else "",
      "\n  index error:", x$obs.error, if(x$catch.error) sprintf("\n  catch error: sdc = %s", if(length(x$sdc) > 1) paste(range(signif(x$sdc, 2)), collapse = "-") else sprintf("%.2f", x$sdc)) else "",
      if(!is.null(x$ebsb)) paste0("\n  Selex: ", length(x$ebsb), " EB/SB curve(s), removals ", if(x$ebsb.removal) "on EB" else "on SB") else "", "\n")
  invisible(x)
}
