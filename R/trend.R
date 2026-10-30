# Data exploration (JARA): state-space trend models for joint indices, census abundance and catch reconstruction, in RTMB
#
# relative: log state x_t = x_t-1 + r + e_t, e_t ~ N(0, sdp^2); log I_it = log q_i + x_t + u_it,
#           u_it ~ N(0, se_it^2 + fixed.obsE^2 + sdi_i^2); q_1 = 1 (the state is on the scale of the first index)
# census:   one state per area (absolute abundance, q = 1), x_kt = x_kt-1 + r_k + e_kt; total N_t = sum_k exp(x_kt)
# timeblock: the drift r differs before and after the block year (JARA timeblock)

.trend_nll <- function(par, dat){
  RTMB::getAll(par, dat)
  dnorm <- RTMB::dnorm
  y <- RTMB::OBS(y)
  nll <- 0
  sdp <- exp(logsdp)
  if(census){
    X <- x; dim(X) <- c(nA, N)
    for(k in seq_len(nA)) nll <- nll - dnorm(X[k, 1], x1.pr[k], 5, log = TRUE)     # diffuse prior on the first state
    for(k in seq_len(nA)) for(t in 2:N) nll <- nll - dnorm(X[k, t], X[k, t - 1] + r[(k - 1)*nb + rb[t]]*drift, sdp, log = TRUE)
    mu <- X[cbind(oi, ot)]
    eX <- exp(X); Ntot <- eX[1, ]; if(nA > 1) for(k in 2:nA) Ntot <- Ntot + eX[k, ]
    logN <- log(Ntot)
    RTMB::ADREPORT(logN); RTMB::REPORT(Ntot)
  } else {
    if(iid){                                       # Conn (2010): independent year effects n_t ~ N(mean, 1)
      nll <- nll - sum(dnorm(x, x1.pr[1], 1, log = TRUE))
    } else {                                       # JARA: random walk with drift
      nll <- nll - dnorm(x[1], x1.pr[1], 5, log = TRUE)                                # diffuse prior on the first state
      for(t in 2:N) nll <- nll - dnorm(x[t], x[t - 1] + r[rb[t]]*drift, sdp, log = TRUE)
    }
    mu <- logq[oi] + x[ot]
    lq <- logq; RTMB::REPORT(lq)
  }
  sdi <- exp(logsdi)
  sdo <- sqrt(se2 + fe2 + (if(vw == 1) sdi[oi] else if(vw == 2) sdi[1] + 0*oi else 0*oi)^2)
  nll <- nll - sum(dnorm(y, mu, sdo, log = TRUE))
  # weak priors: process SD and additional observation SD (lognormal median, log.sd)
  if(!iid) nll <- nll - dnorm(logsdp, log(sdp.pr[1]), sdp.pr[2], log = TRUE)
  if(vw > 0) nll <- nll - sum(dnorm(logsdi, log(sdi.pr[1]), sdi.pr[2], log = TRUE))
  Imean <- mu; sdI <- sdo
  RTMB::ADREPORT(x); RTMB::ADREPORT(Imean)
  RTMB::REPORT(x); RTMB::REPORT(Imean); RTMB::REPORT(sdI); RTMB::REPORT(sdp); RTMB::REPORT(sdi)
  nll
}

#' State-space trend models for abundance indices (JARA)
#'
#' Fits the JARA models in RTMB: a random walk with drift on log abundance observed through one or several indices. The
#' relative joint index follows the hierarchical analysis of multiple noisy abundance indices of Conn (2010), as used in
#' SEDAR assessments; JARA's random-walk form gives practically the same joint index, and \code{state = "iid"} fits
#' Conn's independent-year-effects form directly.
#' \describe{
#'   \item{\code{"relative"}}{joint index of relative abundance indices (CPUE, surveys): each index has its own
#'     catchability q (q of the first index fixed at 1) and observation variance
#'     \eqn{se^2 + fixed.obsE^2 + sdi^2}; the common state is the joint trend.}
#'   \item{\code{"census"}}{indices are absolute abundance (e.g. swept-area or acoustic biomass) of different areas:
#'     one state per area (q = 1), smoothed through missing years, summed to a total abundance with uncertainty.}
#' }
#' \code{variance.weighting}: \code{"model"} estimates an additional observation SD per index, \code{"equal"} one common
#' additional SD, \code{"fixed"} uses only the input SEs plus \code{fixed.obsE}. Weak lognormal priors are placed on the
#' process SD (\code{sdp.pr}) and the additional observation SD (\code{sdi.pr}). The fit carries ordinary and OSA residuals
#' (\code{RTMB::oneStepPredict}, exact for this Gaussian model) in the same format as a \code{jabba3} fit, so
#' \code{jbplot_fits()}, \code{jbplot_residuals()}, \code{jbplot_runs()} and \code{jbplot_ppc()} apply directly. As in
#' JABBA3 the residual plots default to OSA residuals; \code{type = "ordinary"} gives the log residuals used by JARA.
#' @param I data.frame(year, index1, ...) or FLIndices / FLIndexBiomass
#' @param se optional data.frame of log-scale SEs with the same columns (ignored for FLIndices, which carry
#'   \code{index.var})
#' @param type "relative" or "census"
#' @param fixed.obsE additional fixed observation SD added to all indices (as in JARA)
#' @param variance.weighting "model", "equal" or "fixed"
#' @param state relative model only: \code{"rw"} (default), the JARA random walk with drift on the log joint state; or
#'   \code{"iid"}, independent year effects \eqn{n_t \sim N(\bar n, 1)} without time-series structure, i.e. the
#'   hierarchical model of Conn (2010) used in SEDAR assessments (\eqn{\log U_{it} \sim N(n_t + c_i, \sigma_{p,i}^2 +
#'   \sigma_{s,it}^2)}, with \code{variance.weighting = "model"} for the index-specific \eqn{\sigma_{p,i}})
#' @param drift estimate the drift r (random walk with drift); FALSE for a plain random walk
#' @param timeblock optional year(s): the drift (mean rate of change) differs before and after each block year (as JARA's
#'   \code{timeblock}); missing values are interpolated with the block-specific trend
#' @param na.below optional threshold: values below it are treated as unobserved (NA), e.g. incomplete low reports
#' @param sdp.pr,sdi.pr lognormal priors c(median, log.sd) on the process SD and the additional observation SD
#' @param se.na SE used where an index value has no SE
#' @param engine "rtmb" (mode, Laplace) or "nuts" (tmbstan)
#' @param nuts list of tmbstan settings (chains, iter, warmup, seed, cores)
#' @param osa compute OSA residuals
#' @param run label
#' @return object of class \code{jb3trend}
#' @references Conn, P.B. (2010). Hierarchical analysis of multiple noisy abundance indices. Canadian Journal of Fisheries
#'   and Aquatic Sciences, 67, 108-120. doi:10.1139/F09-175
#' @examples
#' \dontrun{
#' data(iccat)
#' tr <- jb3_trend(iccat$bum$cpue, iccat$bum$se, fixed.obsE = 0.1)
#' jbplot_jointindex(tr); jbplot_fits(tr); jbplot_runs(tr, type = "ordinary"); jbplot_ppc(tr)
#' ji <- jb3_jointindex(tr)
#' }
#' @export
jb3_trend <- function(I, se = NULL, type = c("relative", "census"), state = c("rw", "iid"), fixed.obsE = 0,
                      variance.weighting = c("model", "fixed", "equal"), drift = TRUE, timeblock = NULL, na.below = NULL,
                      sdp.pr = c(0.1, 1),
                      sdi.pr = c(0.1, 1), se.na = 0.2, engine = c("rtmb", "nuts"), nuts = list(), osa = TRUE, run = NULL){
  type <- match.arg(type); state <- match.arg(state); variance.weighting <- match.arg(variance.weighting); engine <- match.arg(engine)
  iid <- type == "relative" && state == "iid"
  if(iid){ drift <- FALSE; timeblock <- NULL }
  ix <- .index_df(I, se); Id <- ix$I; Sd <- ix$SE
  has <- Id$year[rowSums(as.matrix(Id[, -1, drop = FALSE]) > 0, na.rm = TRUE) > 0]   # trim leading/trailing years without data
  yrs <- seq(min(has), max(has)); N <- length(yrs)
  Id <- Id[match(yrs, Id$year), ]; Sd <- Sd[match(yrs, Sd$year), ]; Id$year <- Sd$year <- yrs
  nm <- names(Id)[-1]; nI <- length(nm)
  if(!is.null(na.below)) for(i in seq_len(nI)) Id[[i + 1]][is.finite(Id[[i + 1]]) & Id[[i + 1]] < na.below] <- NA
  long <- do.call(rbind, lapply(seq_len(nI), function(i){
    v <- Id[[i + 1]]; s <- Sd[[i + 1]]; ok <- which(is.finite(v) & v > 0)
    if(!length(ok)) return(NULL)
    data.frame(index = nm[i], i = i, t = ok, year = yrs[ok], obs = v[ok], se = ifelse(is.finite(s[ok]), s[ok], se.na))
  }))
  if(is.null(long) || !nrow(long)) stop("no positive index values", call. = FALSE)
  census <- type == "census"
  vw <- c(fixed = 0, model = 1, equal = 2)[[variance.weighting]]
  dat <- list(y = log(long$obs), oi = long$i, ot = long$t, se2 = long$se^2, fe2 = fixed.obsE^2, N = N, nA = nI,
              census = census, vw = vw, drift = as.numeric(drift), sdp.pr = sdp.pr, sdi.pr = sdi.pr, iid = iid)
  tb <- sort(unique(timeblock[timeblock > min(yrs) & timeblock < max(yrs)])); nb <- length(tb) + 1
  dat$nb <- nb; dat$rb <- as.integer(findInterval(yrs, tb, left.open = TRUE) + 1)   # change into year t: block 1 up to the block year, block 2 after
  nd <- if(census) nI else 1
  # initial values: index means; states from the (scaled) mean of the available indices
  lm0 <- tapply(log(long$obs), long$i, mean)
  if(census){
    X0 <- sapply(seq_len(N), function(t) vapply(seq_len(nI), function(k){
      w <- long$i == k; stats::approx(long$t[w], log(long$obs[w]), xout = t, rule = 2)$y }, 0))
    X0 <- matrix(X0, nI, N)
    par <- list(x = as.vector(X0), r = rep(0, nI*nb), logsdp = log(0.1), logsdi = rep(log(0.1), nI))
  } else {
    z <- log(long$obs) - (lm0 - lm0[1])[as.character(long$i)]
    x0 <- stats::approx(long$t, z, xout = seq_len(N), rule = 2, ties = mean)$y
    par <- list(x = x0, logq = as.numeric(lm0 - lm0[1]), r = rep(0, nb), logsdp = log(0.1), logsdi = rep(log(0.1), nI))
  }
  dat$x1.pr <- if(census) X0[, 1] else if(iid) mean(x0) else x0[1]
  map <- list()
  if(!drift) map$r <- factor(rep(NA, length(par$r)))
  if(vw == 0) map$logsdi <- factor(rep(NA, nI)) else if(vw == 2) map$logsdi <- factor(c(1, rep(NA, nI - 1)))
  if(!census) map$logq <- factor(c(NA, seq_len(nI - 1)))
  if(iid) map$logsdp <- factor(NA)
  f <- function(p) .trend_nll(p, dat)
  obj <- RTMB::MakeADFun(f, par, random = "x", map = map, silent = TRUE)
  opt <- suppressWarnings(stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 2000, eval.max = 4000)))
  for(j in 1:2){ g <- as.numeric(obj$gr(opt$par)); h <- tryCatch(stats::optimHess(opt$par, obj$fn, obj$gr), error = function(e) NULL)
    if(is.null(h) || max(abs(g)) < 1e-8) break
    st <- tryCatch(solve(h, g), error = function(e) NULL); if(is.null(st)) break
    np <- opt$par - st; fn <- obj$fn(np); if(is.finite(fn) && fn <= opt$objective) { opt$par <- np; opt$objective <- fn } }
  obj$fn(opt$par)
  sdr <- RTMB::sdreport(obj, opt$par)
  conv <- .conv_tier(obj, opt, sdr)
  rep <- obj$report(obj$env$last.par.best)
  ss <- summary(sdr, "report")
  pick <- function(n) ss[rownames(ss) == n, , drop = FALSE]
  xs <- pick("x"); im <- pick("Imean")
  # states
  if(census){
    st <- data.frame(area = rep(nm, N), year = rep(yrs, each = nI), est = exp(xs[, 1]), se.log = xs[, 2])
    tot <- pick("logN"); total <- data.frame(year = yrs, est = exp(tot[, 1]), se.log = tot[, 2])
    total$lci <- exp(tot[, 1] - 1.96*tot[, 2]); total$uci <- exp(tot[, 1] + 1.96*tot[, 2])
  } else {
    st <- data.frame(area = "trend", year = yrs, est = exp(xs[, 1]), se.log = xs[, 2]); total <- NULL
  }
  st$lci <- st$est*exp(-1.96*st$se.log); st$uci <- st$est*exp(1.96*st$se.log)
  # index fits in the jabba3 format
  sd <- rep$sdI; hat <- exp(im[, 1]); sh <- im[, 2]
  index <- data.frame(index = long$index, year = long$year, obs = long$obs, hat = hat, lci = hat*exp(-1.96*sh), uci = hat*exp(1.96*sh),
                      sd = sd, ppd_lci = hat*exp(-1.96*sqrt(sh^2 + sd^2)), ppd_uci = hat*exp(1.96*sqrt(sh^2 + sd^2)),
                      residual = log(long$obs) - im[, 1], std.residual = (log(long$obs) - im[, 1])/sd, osa = NA, se = long$se)
  if(osa) index$osa <- tryCatch(RTMB::oneStepPredict(obj, method = "fullGaussian", discrete = FALSE, trace = FALSE)$residual,
                                error = function(e) NA)
  fx <- summary(sdr, "fixed")
  pars <- data.frame(par = rownames(fx), est = fx[, 1], se = fx[, 2], row.names = NULL)
  out <- structure(list(type = type, run = run %||% paste("Trend", type), years = yrs, indices = nm, index = index,
                        states = st, total = total, q = if(!census) exp(rep$lq) else rep(1, nI), sdp = rep$sdp,
                        sdi = if(vw == 0) rep(0, nI) else if(vw == 2) rep(rep$sdi[1], nI) else rep$sdi,
                        r = if(drift) matrix(fx[rownames(fx) == "r", 1], nb, nd, dimnames = list(.tb_lab(yrs, tb), if(census) nm else "trend")) else 0,
                        timeblock = tb, pars = pars, convergence = conv, obj = obj, opt = opt,
                        sdr = sdr, data = list(I = Id, SE = Sd, long = long), settings = list(fixed.obsE = fixed.obsE,
                        variance.weighting = variance.weighting, drift = drift, timeblock = tb, na.below = na.below, state = state),
                        engine = engine), class = "jb3trend")
  if(engine == "nuts") out <- .trend_nuts(out, nuts)
  out
}

.tb_lab <- function(yrs, tb){ b <- c(min(yrs), tb, max(yrs)); paste0(c(b[1], tb + 1), "-", b[-1]) }

# draws of a trend fit: MVN at the mode (or NUTS samples), with the index means and SDs for posterior predictive checks
.trend_draws <- function(fit, ndraws = 1000, seed = 1){
  if(!is.null(fit$draws)) return(fit$draws)
  set.seed(seed); obj <- fit$obj
  sdr <- RTMB::sdreport(obj, fit$opt$par, getJointPrecision = TRUE)
  Q <- sdr$jointPrecision; mode <- obj$env$last.par.best
  S <- tryCatch(solve(as.matrix(Q)), error = function(e) NULL)
  if(is.null(S)) stop("joint precision not invertible", call. = FALSE)
  e <- eigen((S + t(S))/2, symmetric = TRUE)
  d <- t(mode + e$vectors %*% (sqrt(pmax(e$values, 0))*matrix(stats::rnorm(length(mode)*ndraws), length(mode))))
  .trend_derive(fit, d, "mvn")
}

.trend_derive <- function(fit, d, type){
  reps <- lapply(seq_len(nrow(d)), function(i) fit$obj$report(d[i, ]))
  list(Ihat = do.call(rbind, lapply(reps, `[[`, "Imean")), sdI = do.call(rbind, lapply(reps, `[[`, "sdI")),
       x = do.call(rbind, lapply(reps, function(r) as.numeric(r$x))), type = type)
}

.trend_nuts <- function(fit, nuts){
  .need("tmbstan", "engine = 'nuts' needs tmbstan")
  st <- utils::modifyList(list(chains = 3, iter = 2000, seed = 1, cores = 1, control = list(adapt_delta = 0.95)), nuts)
  st$warmup <- st$warmup %||% floor(st$iter/2)
  sf <- tmbstan::tmbstan(fit$obj, chains = st$chains, iter = st$iter, warmup = st$warmup, seed = st$seed, cores = st$cores,
                         control = st$control, laplace = FALSE, silent = TRUE, refresh = 0,
                         init = lapply(seq_len(st$chains), function(i) fit$obj$env$last.par.best))
  pm <- as.matrix(sf); pm <- pm[, colnames(pm) != "lp__", drop = FALSE]
  fit$draws <- .trend_derive(fit, pm, "nuts"); fit$nuts <- list(stanfit = sf, settings = st)
  # posterior medians and quantiles of the states
  x <- fit$draws$x; q <- apply(x, 2, stats::quantile, c(0.025, 0.5, 0.975))
  fit$states$est <- exp(q[2, ]); fit$states$lci <- exp(q[1, ]); fit$states$uci <- exp(q[3, ])
  fit$states$se.log <- apply(x, 2, stats::sd)
  fit
}

#' Joint index from a trend fit
#'
#' As JARA's \code{jointindex()}: for each year the precision-weighted mean of the q-scaled observed indices
#' (\eqn{\log I_{it} - \log q_i}, weights 1/total observation variance) with its log-scale SE, and the smoothed state-space
#' trend. Years without observations are NA, as in JARA, or take the state estimate and its SE with \code{fill = TRUE}
#' (column \code{source}). The joint index is scaled
#' to a mean of 1 and can be passed to \code{\link{jb3_data}()} as a single index with \code{log.se}. For a census fit the
#' table holds the total abundance (sum of the area states).
#' Either format goes straight into \code{\link{jb3_data}()}: the data.frame (class \code{jb3jointindex}) as
#' \code{indices} (its \code{log.se} is used as the SE), or the \code{FLIndexBiomass} with \code{index.var} = log.se.
#' @param fit a \code{jb3trend} fit
#' @param format "data.frame" or "FLIndexBiomass"
#' @param name index name
#' @param fill fill years without observations with the state estimate (default FALSE: NA, as JARA)
#' @return data.frame(year, trend, index, log.se, lci, uci, source) of class \code{jb3jointindex}, or an FLIndexBiomass
#' @examples
#' \dontrun{
#' ji <- jb3_jointindex(tr)
#' d  <- jb3_data(catch, ji)                                   # joint index as the single index, log.se as SE
#' d  <- jb3_data(catch, FLCore::FLIndices(joint = jb3_jointindex(tr, "FLIndexBiomass")))
#' }
#' @export
jb3_jointindex <- function(fit, format = c("data.frame", "FLIndexBiomass"), name = "joint", fill = FALSE){
  format <- match.arg(format)
  stopifnot(inherits(fit, "jb3trend"))
  yrs <- fit$years
  if(fit$type == "census"){
    tt <- fit$total
    out <- data.frame(year = yrs, trend = tt$est, index = tt$est, log.se = tt$se.log, lci = tt$lci, uci = tt$uci, source = "state")
  } else {
    l <- fit$data$long; lq <- log(fit$q)[l$i]
    v <- l$se^2 + fit$settings$fixed.obsE^2 + fit$sdi[l$i]^2
    z <- log(l$obs) - lq; w <- 1/v
    agg <- do.call(rbind, lapply(split(data.frame(t = l$t, z = z, w = w), l$t), function(a)
      data.frame(t = a$t[1], m = sum(a$w*a$z)/sum(a$w), s = sqrt(1/sum(a$w)))))
    st <- fit$states
    m <- log(st$est); s <- st$se.log; src <- rep("state", length(yrs))
    m[agg$t] <- agg$m; s[agg$t] <- agg$s; src[agg$t] <- "indices"
    sc <- mean(st$est)
    out <- data.frame(year = yrs, trend = st$est/sc, index = exp(m)/sc, log.se = s, source = src)
    out$lci <- out$index*exp(-1.96*s); out$uci <- out$index*exp(1.96*s)
    out <- out[, c("year", "trend", "index", "log.se", "lci", "uci", "source")]
    if(!fill) out[out$source == "state", c("index", "log.se", "lci", "uci")] <- NA
  }
  if(format == "data.frame") return(structure(out, class = c("jb3jointindex", "data.frame")))
  .need("FLCore", "format = 'FLIndexBiomass' needs FLCore")
  q <- FLCore::FLQuant(out$index, dimnames = list(age = "all", year = yrs))
  fi <- FLCore::FLIndexBiomass(index = q)
  FLCore::index.var(fi) <- FLCore::FLQuant(out$log.se, dimnames = list(age = "all", year = yrs))
  FLCore::name(fi) <- name; fi
}

#' Reconstruct catch series through missing years
#'
#' Census-type state-space smoothing of one or several catch series (e.g. landings by country or fleet), as JARA's census
#' catch reconstruction: a random walk (with optional drift and time blocks) on log catch per series, observed with a
#' near-zero error where the catch is known, so known values are reproduced and missing years are filled with the
#' series-specific (and block-specific) trend; the series are summed to a total with uncertainty.
#'
#' The output goes straight into \code{\link{jb3_data}()} as \code{catch}: the total catch (observed where all series are
#' known, filled otherwise) with its year-specific log-scale SE (\code{se.log}: \code{se} in fully observed years, the SE
#' of the estimated total where any series was filled), used by the catch error of the model in place of the constant
#' \code{sdc}. \code{jb3_fill_table()} gives the original series with the NAs replaced by predictions and the total.
#' @param catch data.frame(year, catch[, more series]); NA for missing values (years can also be absent)
#' @param se observation SD (log scale) for known catches (JARA: \code{fixed.obsE = 0.001})
#' @param drift random walk with drift (default TRUE with a time block, otherwise FALSE)
#' @param timeblock optional year(s) after which the drift changes (JARA \code{timeblock})
#' @param na.below optional threshold: catches below it are treated as unobserved
#' @param years optional full year range (e.g. to extend backwards)
#' @param sdp.pr lognormal prior on the process SD
#' @return data.frame(year, obs (sum of observed series, NA if any missing), catch, se.log, est, lci, uci, filled) of class
#'   \code{jb3catchfill}, with the census fit in attribute \code{"fit"}
#' @examples
#' \dontrun{
#' cc <- iccat$bet$catch; cc[c(5:8, 30:31), 2] <- NA
#' cf <- jb3_catchfill(cc)
#' d <- jb3_data(cf, iccat$bet$cpue, iccat$bet$se)   # filled catch with year-specific catch SE
#' fit <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "fox"))
#' }
#' @export
jb3_catchfill <- function(catch, se = 0.01, drift = !is.null(timeblock), timeblock = NULL, na.below = NULL, years = NULL,
                          sdp.pr = c(0.2, 1)){
  cc <- as.data.frame(catch); cc[[1]] <- as.numeric(cc[[1]]); names(cc)[1] <- "year"
  if(is.null(names(cc)) || ncol(cc) == 2) names(cc)[2] <- names(cc)[2] %||% "catch"
  yrs <- years %||% seq(min(cc$year), max(cc$year))
  I <- data.frame(year = yrs); S <- data.frame(year = yrs)
  for(v in names(cc)[-1]){ I[[v]] <- as.numeric(cc[[v]])[match(yrs, cc$year)]; S[[v]] <- se }
  if(!is.null(na.below)) for(v in names(I)[-1]) I[[v]][is.finite(I[[v]]) & I[[v]] < na.below] <- NA
  fit <- jb3_trend(I, S, type = "census", variance.weighting = "fixed", drift = drift, timeblock = timeblock, sdp.pr = sdp.pr,
                   osa = FALSE, run = "Catch reconstruction")
  ft <- jb3_fill_table(fit)
  obs <- I[, -1, drop = FALSE]; anymiss <- rowSums(!is.finite(as.matrix(obs)) | as.matrix(obs) <= 0) > 0
  tt <- fit$total
  out <- data.frame(year = yrs, obs = ifelse(anymiss, NA, rowSums(as.matrix(obs))), catch = ft$Total,
                    se.log = ifelse(anymiss, tt$se.log, se), est = tt$est, lci = tt$lci, uci = tt$uci, filled = anymiss)
  structure(out, fit = fit, class = c("jb3catchfill", "data.frame"))
}

#' Original series with missing values filled
#'
#' For a census trend fit (or \code{jb3_catchfill()} output): the original series with NAs (and values treated as
#' unobserved) replaced by the fitted states, and their total (as JARA's reconstructed catch table).
#' @param x a census \code{jb3trend} fit or \code{jb3catchfill} output
#' @param digits rounding
#' @return data.frame(year, series..., Total)
#' @export
jb3_fill_table <- function(x, digits = NULL){
  fit <- if(inherits(x, "jb3catchfill")) attr(x, "fit") else x
  if(!inherits(fit, "jb3trend") || fit$type != "census") stop("needs a census trend fit or jb3_catchfill() output", call. = FALSE)
  I <- fit$data$I; st <- fit$states; out <- data.frame(year = fit$years)
  for(v in fit$indices){ o <- I[[v]]; p <- st$est[st$area == v]
    out[[v]] <- ifelse(is.finite(o) & o > 0, o, p) }
  out$Total <- rowSums(as.matrix(out[, -1, drop = FALSE]))
  if(!is.null(digits)) out[, -1] <- round(out[, -1], digits)
  out
}

#' @export
print.jb3trend <- function(x, ...){
  cat("JABBA3 trend fit (", x$type, "): ", x$run, "\n", sep = "")
  cat("  indices:", paste(x$indices, collapse = ", "), "| years", min(x$years), "-", max(x$years), "\n")
  cat("  convergence:", x$convergence$tier, sprintf("(max gradient %.1e)", x$convergence$maxgrad), "\n")
  if(identical(x$settings$state, "iid")) cat("  state: independent year effects (Conn 2010)\n") else cat(sprintf("  process SD = %.3f", x$sdp), if(length(x$timeblock)) paste("| time block after", paste(x$timeblock, collapse = ", ")), "\n  drift r:\n")
  if(!identical(x$settings$state, "iid")) print(round(x$r, 3))
  if(x$type == "relative") print(data.frame(index = x$indices, q = signif(x$q, 3), sdi = signif(x$sdi, 3)), row.names = FALSE)
  invisible(x)
}

#' Trend plots
#'
#' \code{jbplot_trend()}: relative fits, the q-scaled observed indices (by index) with the state-space trend and its 95\%
#' interval; census fits, the area states with the total abundance. \code{jbplot_jointindex()}: the joint index with 95\%
#' intervals and the state-space trend (red dashed), as JARA's \code{jrplot_jointindex()}. \code{jbplot_catchfill()}:
#' observed and reconstructed catch.
#' @param fit a \code{jb3trend} fit (or \code{jb3_catchfill()} output for \code{jbplot_catchfill})
#' @return a ggplot
#' @export
jbplot_trend <- function(fit){
  stopifnot(inherits(fit, "jb3trend"))
  st <- fit$states; l <- fit$data$long
  if(fit$type == "census"){
    tt <- fit$total; tt$area <- "Total"
    d <- rbind(st[, c("area", "year", "est", "lci", "uci")], tt[, c("area", "year", "est", "lci", "uci")])
    d$area <- factor(d$area, c(fit$indices, "Total"))
    o <- data.frame(area = factor(l$index, levels(d$area)), year = l$year, obs = l$obs)
    ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$est)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey75", alpha = 0.6) +
      ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(data = o, ggplot2::aes(.data$year, .data$obs), shape = 21, fill = "white") +
      ggplot2::facet_wrap(~area, scales = "free_y") + ggplot2::expand_limits(y = 0) + scale_x_year() +
      ggplot2::labs(x = "Year", y = "Abundance", subtitle = "Census: area states and total (95% intervals)") + theme_jabba3()
  } else {
    sc <- mean(st$est); o <- data.frame(index = factor(l$index, fit$indices), year = l$year, obs = l$obs/fit$q[l$i]/sc)
    ggplot2::ggplot(st, ggplot2::aes(.data$year, .data$est/sc)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci/sc, ymax = .data$uci/sc), fill = "grey75", alpha = 0.6) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(data = o, ggplot2::aes(.data$year, .data$obs, fill = .data$index), shape = 21, size = 1.8, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = jb3cols(max(2, length(fit$indices)))[seq_along(fit$indices)], name = NULL) +
      ggplot2::expand_limits(y = 0) + scale_x_year() +
      ggplot2::labs(x = "Year", y = "Relative abundance", subtitle = "Indices scaled by q; state-space trend with 95% interval") +
      theme_jabba3() + ggplot2::theme(legend.position = "right")
  }
}

#' @rdname jbplot_trend
#' @export
jbplot_jointindex <- function(fit){
  ji <- as.data.frame(jb3_jointindex(fit))
  ggplot2::ggplot(ji, ggplot2::aes(.data$year)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$trend), colour = "red", linetype = 2, linewidth = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), width = 0.3, colour = "grey30") +
    ggplot2::geom_line(ggplot2::aes(y = .data$index), linewidth = 0.2) +
    ggplot2::geom_point(ggplot2::aes(y = .data$index, shape = .data$source), fill = "white", size = 2) +
    ggplot2::scale_shape_manual(values = c(indices = 21, state = 24), name = NULL) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = if(fit$type == "census") "Total abundance" else "Joint index",
                  subtitle = "Joint index with 95% intervals; red dashed: state-space trend") + theme_jabba3()
}

#' @rdname jbplot_trend
#' @export
jbplot_catchfill <- function(fit){
  d <- as.data.frame(fit)
  ggplot2::ggplot(d, ggplot2::aes(.data$year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = .data$catch), linewidth = 0.6) +
    ggplot2::geom_point(data = d[!d$filled, ], ggplot2::aes(y = .data$obs), size = 1.4) +
    ggplot2::geom_point(data = d[d$filled, ], ggplot2::aes(y = .data$catch), shape = 21, fill = "orange", size = 2.2) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Catch", subtitle = "Observed (black) and reconstructed (orange) catch, 95% interval") + theme_jabba3()
}
