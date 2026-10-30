# Draws: MVN pseudo-posterior at the mode, and NUTS via tmbstan

# derived quantities from full parameter vectors (rows = draws)
.derive <- function(fit, draws){
  obj <- fit$obj; mu <- fit$mucatch; yrs <- fit$years
  reps <- lapply(seq_len(nrow(draws)), function(i) obj$report(draws[i, ]))
  col <- function(nm, sc = 1){ x <- do.call(rbind, lapply(reps, function(r) as.numeric(r[[nm]])*sc)); colnames(x) <- yrs; x }
  sc <- function(nm, sc = 1) vapply(reps, function(r) as.numeric(r[[nm]])[1]*sc, 0)
  K <- sc("K", mu)
  out <- list(B = col("P")*K, `B/Bmsy` = col("BBmsy"), F = col("H"), `F/Fmsy` = col("HHmsy"), `B/K` = col("P"),
              Catch = col("estC", mu), `Process deviation` = col("procdev"),
              pars = data.frame(K = K, r = sc("r"), bkinit = sc("bk"), sdb = sc("sdb"), m = sc("m"), kdyn = sc("kdyn"),
                                Bmsy = sc("BmsyK")*K, Fmsy = sc("Hmsy"), MSY = sc("MSY", mu)),
              q = do.call(rbind, lapply(reps, function(r) r$q/mu)), sdi = do.call(rbind, lapply(reps, function(r) r$sdi)),
              Ihat = do.call(rbind, lapply(reps, function(r) r$Imean)), sdI = do.call(rbind, lapply(reps, function(r) r$sdI)),
              par = draws)
  out
}

#' Pseudo-posterior draws at the mode
#'
#' Draws all fixed and random effects from the multivariate normal with the joint precision at the mode (Laplace
#' approximation), and returns derived trajectories and parameters.
#' @param fit a \code{jabba3} fit (rtmb)
#' @param ndraws number of draws
#' @param seed random seed
#' @export
jb3_draws <- function(fit, ndraws = 1000, seed = 1){
  set.seed(seed)
  obj <- fit$obj
  sdr <- RTMB::sdreport(obj, fit$opt$par, getJointPrecision = TRUE)
  Q <- sdr$jointPrecision; mode <- obj$env$last.par.best
  z <- matrix(stats::rnorm(nrow(Q)*ndraws), nrow(Q))
  L <- tryCatch(Matrix::Cholesky(methods::as(Q, "CsparseMatrix"), super = TRUE, perm = TRUE), error = function(e) NULL)
  x <- if(!is.null(L)) as.matrix(Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "Pt")) else {
    S <- solve(as.matrix(Q)); e <- eigen((S + t(S))/2, symmetric = TRUE); e$vectors %*% (sqrt(pmax(e$values, 0))*z) }
  d <- t(mode + x)
  out <- .derive(fit, d); out$type <- "mvn"
  out
}

# NUTS with tmbstan on the same objective (all random effects sampled)
jb3_nuts <- function(fit, nuts = list()){
  .need("tmbstan", "engine = 'nuts' needs tmbstan")
  st <- utils::modifyList(list(chains = 3, iter = 4000, seed = 1, cores = 1, thin = 1,
                               control = list(adapt_delta = 0.95, max_treedepth = 12)), nuts)
  if(is.null(st$warmup)) st$warmup <- floor(st$iter/2)               # as rstan: half of iter
  if(st$warmup >= st$iter) stop("nuts: warmup must be smaller than iter", call. = FALSE)
  obj <- fit$obj
  init <- lapply(seq_len(st$chains), function(i) obj$env$last.par.best)
  sf <- tmbstan::tmbstan(obj, chains = st$chains, iter = st$iter, warmup = st$warmup, seed = st$seed, cores = st$cores,
                         thin = st$thin, control = st$control, laplace = FALSE, init = init, silent = TRUE, refresh = 0)
  pm <- as.matrix(sf); pm <- pm[, colnames(pm) != "lp__", drop = FALSE]
  fit$draws <- .derive(fit, pm); fit$draws$type <- "nuts"
  sm <- rstan::summary(sf)$summary
  fit$nuts <- list(stanfit = sf, rhat.max = max(sm[, "Rhat"], na.rm = TRUE), ess.min = min(sm[, "n_eff"], na.rm = TRUE),
                   divergent = sum(rstan::get_divergent_iterations(sf)), settings = st)
  # convergence check for final runs: rank-normalised R-hat, bulk/tail ESS per chain, divergences
  dg <- tryCatch(jb3_nuts_diags(fit, states = TRUE), error = function(e) NULL)
  if(!is.null(dg)){
    fit$nuts$diags <- dg
    msg <- c(if(max(dg$pars$Rhat, na.rm = TRUE) > 1.01) sprintf("R-hat up to %.3f (> 1.01)", max(dg$pars$Rhat, na.rm = TRUE)),
             if(min(dg$pars$ess_bulk, na.rm = TRUE) < 100*st$chains) sprintf("bulk ESS down to %.0f (< 100 per chain)", min(dg$pars$ess_bulk, na.rm = TRUE)),
             if(min(dg$pars$ess_tail, na.rm = TRUE) < 100*st$chains) sprintf("tail ESS down to %.0f (< 100 per chain)", min(dg$pars$ess_tail, na.rm = TRUE)),
             if(sum(dg$sampler$divergent) > 0) sprintf("%d divergent transitions", sum(dg$sampler$divergent)),
             if(any(dg$sampler$EBFMI < 0.3)) "E-BFMI < 0.3")
    if(length(msg)) warning("NUTS not converged: ", paste(msg, collapse = "; "),
                            ". Run longer (nuts = list(iter = ...)) and check jbplot_nutsdiags().", call. = FALSE)
  }
  # posterior medians and quantiles replace the mode summaries
  fit$ts <- .ts_from_draws(fit, fit$draws)
  pp <- fit$draws$pars
  qq <- function(x) c(stats::median(x), stats::quantile(x, c(0.025, 0.975), names = FALSE))
  for(nm in intersect(c("K", "r", "bkinit", "sdb", "m", "kdyn"), rownames(fit$pars))) fit$pars[nm, c("est", "lci", "uci")] <- qq(pp[[nm]])
  for(nm in c("Bmsy", "Fmsy", "MSY", "K")) fit$refpts[nm, c("est", "lci", "uci")] <- qq(pp[[nm]])
  fit
}

.ts_from_draws <- function(fit, dr, probs = c(0.025, 0.5, 0.975)){
  keep <- fit$ts[!fit$ts$quantity %in% c("B", "B/Bmsy", "F", "F/Fmsy", "B/K", "Catch", "Process deviation"), ]
  new <- do.call(rbind, lapply(c("B", "F", "B/Bmsy", "F/Fmsy", "B/K", "Catch", "Process deviation"), function(q){
    x <- dr[[q]]; qs <- apply(x, 2, stats::quantile, probs, na.rm = TRUE)
    data.frame(year = fit$years, quantity = q, est = qs[2, ], lci = qs[1, ], uci = qs[3, ], se.log = NA, run = NA_character_)
  }))
  rbind(new, keep)
}
