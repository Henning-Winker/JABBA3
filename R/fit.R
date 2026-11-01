#' Fit a JABBA3 model
#'
#' Fits the Pella-Tomlinson state-space surplus production model in RTMB (posterior mode of the fixed effects with the
#' log-depletion states Laplace-integrated) or with NUTS (\code{tmbstan}, all states sampled). Optionally runs one of the
#' BioProcDyn workflows (Tier 1 lag, Tier 2 tuned dynamic process error, Selex) automatically.
#'
#' @param data \code{\link{jb3_data}()} object
#' @param priors \code{\link{jb3_priors}()} object (default: JABBA defaults with the FishLife sdb prior)
#' @param model \code{\link{jb3_model}()} object
#' @param engine "rtmb" (default) or "nuts"
#' @param ndraws number of pseudo-posterior draws from the joint precision at the mode (rtmb; 0 = none). Draws feed the
#'   Kobe cloud and projections; intervals come from the delta method otherwise.
#' @param nuts list of tmbstan settings: chains (3), iter (4000), warmup (iter/2), seed, cores, thin, control
#'   (adapt_delta 0.95, max_treedepth 12). NUTS is meant for final runs: convergence is checked after sampling and a warning
#'   is issued if R-hat > 1.01, bulk/tail ESS < 100 per chain, divergences occur or E-BFMI < 0.3 (table in \code{fit$nuts$diags})
#' @param osa compute one-step-ahead residuals (particle filter at the mode)
#' @param workflow "none" (default, manual), "tier1", "tier2", "selex", "selex_tier1", "selex_tier2" (see
#'   \code{\link{jb3_workflow}()})
#' @param amat age at 50\% maturity (production lag for the Tier 1 / Tier 2 workflows)
#' @param run name of the run (scenario)
#' @param se compute standard errors (sdreport). \code{FALSE} gives point estimates only, about twice as fast, for
#'   simulation testing; the convergence verdict then rests on the gradient (no Hessian check)
#' @param peels,index.peels internal: retrospective (all data) and hindcast (index only) peels
#' @param nstart maximum number of optimisation starts
#' @param silent suppress RTMB tracing
#' @return list of class \code{jabba3}
#' @examples
#' data(iccat)
#' d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
#' fit <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"))
#' fit
#' @export
jabba3 <- function(data, priors = jb3_priors(data), model = jb3_model(), engine = c("rtmb", "nuts"), ndraws = 0,
                   nuts = list(), osa = TRUE, workflow = c("none", "tier1", "tier2", "selex", "selex_tier1", "selex_tier2"),
                   amat = NULL, run = NULL, se = TRUE, peels = 0, index.peels = 0, nstart = 4, silent = TRUE, prior.only = FALSE){
  .need("RTMB", "JABBA3 needs RTMB")
  if(!inherits(data, "jb3data")) stop("data must be a jb3_data() object", call. = FALSE)
  if(missing(engine) && (ncol(data$I) == 0 || isTRUE(prior.only))){
    engine <- "nuts"
    message(if(isTRUE(prior.only)) "prior.only = TRUE" else "no abundance index (catch-only)", ": engine = \"nuts\" (set engine = \"rtmb\" for the posterior mode)")
  }
  engine <- match.arg(engine); workflow <- match.arg(workflow)
  if(ncol(data$I) == 0 || isTRUE(prior.only)) osa <- FALSE
  if(workflow != "none") return(jb3_workflow(data, priors, model, workflow = workflow, amat = amat, engine = engine,
                                             ndraws = ndraws, nuts = nuts, osa = osa, se = se))
  ptm <- proc.time()
  setup <- jb3_setup(data, priors, model, peels = peels, index.peels = index.peels, prior.only = prior.only)
  opt <- jb3_optimise(setup, nstart = nstart, silent = silent)
  fit <- jb3_output(opt, setup, sdreport = se || engine == "nuts" || ndraws > 0)
  class(fit) <- c("jabba3", "list")
  fit$run <- run %||% data$stock
  fit$priors <- priors; fit$model <- model; fit$engine <- engine; fit$peels <- peels; fit$index.peels <- index.peels
  fit$prior.only <- isTRUE(prior.only)
  if(engine == "nuts"){
    fit <- jb3_nuts(fit, nuts)
  } else if(ndraws > 0) fit$draws <- jb3_draws(fit, ndraws)
  if(osa && index.peels == 0 && peels == 0) fit$osa <- tryCatch(jb3_osa(fit), error = function(e){ warning("OSA residuals failed: ", conditionMessage(e), call. = FALSE); NULL })
  fit$index <- .index_table(fit)
  fit$auxiliary <- .aux_table(fit)
  fit$time <- (proc.time() - ptm)[["elapsed"]]
  fit
}

#---------------------------------------------------------------------------------------------------------------
# Optimisation as in FLicc: nlminb, Newton polishing, gradient-guarded restarts, Newton-decrement verdict
jb3_optimise <- function(setup, nstart = 4, silent = TRUE){
  dat <- setup$dat
  f <- function(p) jb3_nll(p, dat)
  newton <- function(obj, opt, steps = 3){
    for(j in seq_len(steps)){
      g <- tryCatch(as.numeric(obj$gr(opt$par)), error = function(e) NA)
      if(all(is.finite(g)) && max(abs(g)) < 1e-6) break          # already converged: skip the Hessian
      h <- tryCatch(stats::optimHess(opt$par, obj$fn, obj$gr), error = function(e) NULL)
      if(is.null(h) || any(!is.finite(g))) break
      st <- tryCatch(solve(h, g), error = function(e) NULL); if(is.null(st)) break
      np <- opt$par - st; fn <- tryCatch(obj$fn(np), error = function(e) NA)
      if(is.finite(fn) && fn <= opt$objective + 1e-10){ opt$par <- np; opt$objective <- fn } else break
      if(max(abs(g)) < 1e-8) break
    }
    opt
  }
  run1 <- function(par){
    obj <- tryCatch(RTMB::MakeADFun(f, par, random = setup$random, map = setup$map, silent = silent), error = function(e) NULL)
    if(is.null(obj)) return(NULL)
    opt <- tryCatch(suppressWarnings(stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 2000))),
                    error = function(e) NULL)
    if(is.null(opt) || !is.finite(opt$objective)) return(NULL)
    opt <- newton(obj, opt)
    opt$maxgrad <- max(abs(obj$gr(opt$par)))
    list(obj = obj, opt = opt)
  }
  best <- run1(setup$par); nst <- 1
  ok <- function(b) !is.null(b) && b$opt$convergence == 0 && b$opt$maxgrad < 1e-3
  while(!ok(best) && nst < nstart){
    nst <- nst + 1
    p <- jb3_inits(dat, jitter = 1)
    fix <- setdiff(names(setup$par), c("logK", "logr", "logq", "logP"))
    p[fix] <- setup$par[fix]
    b <- run1(p)
    # guarded: accept only if the objective does not rise and the gradient falls
    if(is.null(best) || (!is.null(b) && b$opt$objective <= best$opt$objective + 1e-6 &&
                         (b$opt$maxgrad < best$opt$maxgrad || b$opt$objective < best$opt$objective - 1e-3))) best <- b
  }
  if(is.null(best)) stop("optimisation failed from all starts", call. = FALSE)
  best$opt <- newton(best$obj, best$opt, 2)
  best$obj$fn(best$opt$par)
  best$nstart <- nst
  best
}

#' Convergence checks
#'
#' FLicc-style convergence table: optimiser code, maximum gradient, positive-definite Hessian, standard errors and the
#' Newton decrement 0.5 g' H^-1 g, with an overall verdict: "ok" (decrement < 0.001), "warning" (<= 0.1) or "fail"
#' (above 0.1, or Hessian not positive definite).
#' @param fit a \code{jabba3} fit
#' @return data.frame of checks with a plain-language interpretation; verdict in attribute "tier"
#' @export
jb3_convergence <- function(fit){
  cv <- fit$convergence
  tab <- data.frame(
    check = c("optimiser code", "max gradient", "Hessian positive definite", "standard errors", "Newton decrement"),
    value = c(cv$code, signif(cv$maxgrad, 3), cv$pdHess, cv$se.ok, signif(cv$decrement, 3)),
    pass = c(cv$code == 0, cv$maxgrad < 1e-3, isTRUE(cv$pdHess), isTRUE(cv$se.ok), isTRUE(cv$decrement < 1e-3)),
    meaning = c("nlminb reported convergence (0)", "largest absolute gradient at the optimum (< 0.001)",
                "the mode is a proper minimum (curvature in every direction)", "all delta-method SEs finite",
                "expected further decrease of the objective by a Newton step (< 0.001 ok, <= 0.1 warning)"))
  attr(tab, "tier") <- cv$tier
  tab
}

.conv_tier <- function(obj, opt, sdr){
  g <- tryCatch(as.numeric(obj$gr(opt$par)), error = function(e) rep(NA, length(opt$par)))
  if(is.null(sdr)){   # point-estimate mode: verdict from the optimiser and the gradient only
    tier <- if(opt$convergence == 0 && all(is.finite(g)) && max(abs(g)) < 1e-3) "ok" else if(all(is.finite(g)) && max(abs(g)) < 0.1) "warning" else "fail"
    return(list(code = opt$convergence, message = opt$message, maxgrad = max(abs(g)), pdHess = NA, se.ok = NA,
                decrement = NA, tier = tier, nll = opt$objective))
  }
  # Newton decrement from the sdreport covariance (inverse Hessian of the marginal), no extra Hessian evaluation
  dec <- if(any(!is.finite(g))) NA else tryCatch(0.5*sum(g*(sdr$cov.fixed %*% g)), error = function(e) NA)
  pd <- isTRUE(sdr$pdHess)
  se.ok <- all(is.finite(suppressWarnings(sqrt(diag(sdr$cov.fixed)))))
  tier <- if(!pd || !is.finite(dec) || dec > 0.1) "fail" else if(dec >= 1e-3) "warning" else "ok"
  list(code = opt$convergence, message = opt$message, maxgrad = max(abs(g)), pdHess = pd, se.ok = se.ok,
       decrement = dec, tier = tier, nll = opt$objective)
}

#---------------------------------------------------------------------------------------------------------------
# Lean output at the mode
jb3_output <- function(best, setup, sdreport = TRUE){
  obj <- best$obj; opt <- best$opt; dat <- setup$dat; mu <- setup$mucatch; yrs <- setup$years; N <- dat$N
  rep <- obj$report(obj$env$last.par.best)
  sdr <- if(sdreport) tryCatch(suppressWarnings(RTMB::sdreport(obj, opt$par)), error = function(e) NULL) else NULL
  conv <- .conv_tier(obj, opt, sdr); conv$nstart <- best$nstart
  g <- function(nm){
    if(is.null(sdr)) return(NULL); i <- names(sdr$value) == nm
    if(!any(i)) return(NULL); cbind(est = sdr$value[i], se = sdr$sd[i])
  }
  lnci <- function(x, scale = 1, est = NULL){
    e <- if(is.null(x)) log(est) else x[, 1]; s <- if(is.null(x)) NA else x[, 2]
    data.frame(est = exp(e)*scale, lci = exp(e - 1.96*s)*scale, uci = exp(e + 1.96*s)*scale, se.log = s)
  }
  # time series (long)
  tsq <- function(q, x, scale = 1, est = NULL) cbind(data.frame(year = yrs, quantity = q), lnci(x, scale, est))
  lp <- if(!is.null(sdr)) sdr$par.random[names(sdr$par.random) == "logP"] else log(rep$P)
  lps <- if(!is.null(sdr)) sqrt(sdr$diag.cov.random[names(sdr$par.random) == "logP"]) else NA
  B <- rep$P*rep$K
  remv <- if(dat$eb.rem) rep$rem else rep(1, N)
  SP <- c((B[-1] - B[-N] + rep$estC[-N]/remv[-N])*remv[-N], NA)*mu
  pd <- g("procdev")
  ts <- rbind(tsq("B", g("logB"), mu, B), tsq("F", g("logF"), 1, rep$H), tsq("B/Bmsy", g("logBBmsy"), 1, rep$BBmsy),
              tsq("F/Fmsy", g("logFFmsy"), 1, rep$HHmsy),
              cbind(data.frame(year = yrs, quantity = "B/K"), data.frame(est = exp(lp), lci = exp(lp - 1.96*lps), uci = exp(lp + 1.96*lps), se.log = lps)),
              data.frame(year = yrs, quantity = "Catch", est = rep$estC*mu,
                         lci = if(dat$catch.err && !is.null(g("estC"))) (g("estC")[, 1] - 1.96*g("estC")[, 2])*mu else rep$estC*mu,
                         uci = if(dat$catch.err && !is.null(g("estC"))) (g("estC")[, 1] + 1.96*g("estC")[, 2])*mu else rep$estC*mu, se.log = NA),
              data.frame(year = yrs, quantity = "Process deviation", est = rep$procdev,
                         lci = if(is.null(pd)) NA else pd[, 1] - 1.96*pd[, 2], uci = if(is.null(pd)) NA else pd[, 1] + 1.96*pd[, 2], se.log = NA),
              data.frame(year = yrs, quantity = "Process SD", est = rep$sdt, lci = NA, uci = NA, se.log = NA),
              data.frame(year = yrs, quantity = "Surplus production", est = SP, lci = NA, uci = NA, se.log = NA))
  if(dat$eb.on){
    EB <- B*remv*mu
    ts <- rbind(ts, data.frame(year = yrs, quantity = "EB", est = EB, lci = NA, uci = NA, se.log = NA))
  }
  ts$run <- NA_character_
  # parameters
  p1 <- function(nm, x, est, scale = 1) cbind(data.frame(par = nm), lnci(x, scale, est))
  q <- g("logq"); sdi <- g("logsdi")
  pars <- rbind(p1("K", g("logK"), rep$K, mu), p1("r", g("logr"), rep$r), p1("bkinit", g("logbk"), rep$bk),
                p1("sdb", if(dat$est.sdb) g("logsdb") else NULL, rep$sdb),
                if(dat$nI > 0) do.call(rbind, lapply(seq_len(dat$nq), function(j) p1(paste0("q", j), if(is.null(q)) NULL else q[j, , drop = FALSE], rep$q[j], 1/mu))),
                if(dat$nI > 0) do.call(rbind, lapply(seq_len(dat$nvar), function(j) p1(paste0("sdi", j), if(is.null(sdi)) NULL else sdi[j, , drop = FALSE], rep$sdi[j]))),
                p1("m", if(dat$est.m) g("logm") else NULL, rep$m))
  if(dat$nA > 0){
    qa <- g("logqA"); sa <- g("logsdA")
    pars <- rbind(pars, do.call(rbind, lapply(seq_len(dat$nA), function(j) rbind(
      p1(paste0("qA.", dat$aux.names[j]), if(is.null(qa) || !dat$est.qA[j]) NULL else qa[j, , drop = FALSE], rep$qA[j]),
      p1(paste0("sdA.", dat$aux.names[j]), if(is.null(sa) || !dat$est.sdA[j]) NULL else sa[j, , drop = FALSE], rep$sdA[j])))))
    if(dat$est.creep){ k <- g("creep"); pars <- rbind(pars, data.frame(par = "creep", est = rep$creep,
      lci = if(is.null(k)) NA else k[1] - 1.96*k[2], uci = if(is.null(k)) NA else k[1] + 1.96*k[2], se.log = NA)) }
  }
  if(dat$dyn){ k <- g("kdyn"); pars <- rbind(pars, data.frame(par = "kdyn", est = rep$kdyn,
    lci = if(is.null(k)) NA else k[1] - 1.96*k[2], uci = if(is.null(k)) NA else k[1] + 1.96*k[2], se.log = NA)) }
  if(dat$est.rho || dat$rho.fix > 0){ k <- g("rho"); pars <- rbind(pars, data.frame(par = "rho", est = rep$rho,
    lci = if(is.null(k)) NA else k[1] - 1.96*k[2], uci = if(is.null(k)) NA else k[1] + 1.96*k[2], se.log = NA)) }
  refpts <- rbind(p1("Bmsy", g("logBmsy"), rep$BmsyK*rep$K, mu), p1("Fmsy", g("logFmsy"), rep$Hmsy), p1("MSY", g("logMSY"), rep$MSY, mu),
                  p1("K", g("logK"), rep$K, mu), p1("BmsyK", NULL, rep$BmsyK))
  rownames(pars) <- pars$par; rownames(refpts) <- refpts$par
  if(!is.null(sdr) && !isTRUE(conv$pdHess)){   # no valid curvature: intervals are not meaningful
    ts[, c("lci", "uci", "se.log")] <- NA; pars[, c("lci", "uci", "se.log")] <- NA; refpts[, c("lci", "uci", "se.log")] <- NA
    warning("Hessian not positive definite: intervals set to NA (convergence: ", conv$tier, ")", call. = FALSE)
  }
  list(stock = setup$data$stock, years = yrs, data = setup$data, ts = ts, pars = pars, refpts = refpts, convergence = conv,
       report = rep, obj = obj, opt = opt, sdr = sdr, setup = setup, mucatch = mu)
}

# index fits table: obs, fitted (with CI), predictive interval, ordinary and OSA residuals
.index_table <- function(fit){
  s <- fit$setup; d <- s$dat; rep <- fit$report
  if(!length(d$obs_t)) return(data.frame(index = character(0), year = numeric(0), obs = numeric(0), hat = numeric(0), lci = numeric(0),
                                         uci = numeric(0), sd = numeric(0), ppd_lci = numeric(0), ppd_uci = numeric(0),
                                         residual = numeric(0), std.residual = numeric(0), osa = numeric(0)))
  im <- if(!is.null(fit$sdr)){ i <- names(fit$sdr$value) == "Imean"; cbind(fit$sdr$value[i], fit$sdr$sd[i]) } else cbind(rep$Imean, NA)
  out <- data.frame(index = s$index.names[d$obs_i], year = s$years[d$obs_t], obs = exp(d$logIobs),
                    hat = exp(im[, 1]), lci = exp(im[, 1] - 1.96*im[, 2]), uci = exp(im[, 1] + 1.96*im[, 2]),
                    sd = rep$sdI, ppd_lci = exp(im[, 1] - 1.96*sqrt(rep$sdI^2 + im[, 2]^2)),
                    ppd_uci = exp(im[, 1] + 1.96*sqrt(rep$sdI^2 + im[, 2]^2)))
  out$residual <- log(out$obs) - log(out$hat)
  out$std.residual <- out$residual/out$sd
  out$osa <- if(!is.null(fit$osa)) fit$osa else NA
  out
}

# auxiliary fits table: obs, fitted (with CI), predictive interval, residuals
.aux_table <- function(fit){
  s <- fit$setup; d <- s$dat; rep <- fit$report
  if(is.null(d$nA) || d$nA == 0) return(NULL)
  am <- if(!is.null(fit$sdr) && any(names(fit$sdr$value) == "Amean")){ i <- names(fit$sdr$value) == "Amean"; cbind(fit$sdr$value[i], fit$sdr$sd[i]) } else cbind(rep$Amean, NA)
  if(!is.null(fit$draws$Ahat)){ am[, 1] <- apply(fit$draws$Ahat, 2, stats::median); am[, 2] <- apply(fit$draws$Ahat, 2, stats::sd) }
  out <- data.frame(series = d$aux.names[d$aobs_i], type = s$data$A.type[d$aobs_i], year = s$years[d$aobs_t], obs = exp(d$logAobs),
                    hat = exp(am[, 1]), lci = exp(am[, 1] - 1.96*am[, 2]), uci = exp(am[, 1] + 1.96*am[, 2]), sd = rep$sdAo,
                    ppd_lci = exp(am[, 1] - 1.96*sqrt(rep$sdAo^2 + am[, 2]^2)), ppd_uci = exp(am[, 1] + 1.96*sqrt(rep$sdAo^2 + am[, 2]^2)))
  out$residual <- log(out$obs) - log(out$hat)
  out$std.residual <- out$residual/out$sd
  out
}
