#' Convert a JABBA input object to JABBA3
#'
#' Translates a \code{JABBA::build_jabba()} object into the three JABBA3 inputs: \code{data}, \code{priors} and
#' \code{model}. Priors are carried over on the JABBA scale (r, K, psi -> bkinit, shape, q1, status prior), as are the
#' lag, catch error, catchability and variance groups, the process-error start, the BioProcDyn options (lognormal
#' sigma.proc prior, fixed sigma.proc, depletion-dependent floor, bias correction) and the Selex EB/SB curves.
#'
#' Differences that cannot be carried over one to one are reported as messages:
#' \itemize{
#'   \item JABBA's inverse-gamma prior on the process variance is replaced by the JABBA3 default (FishLife lognormal
#'     prior on sdb), because the inverse-gamma mode collapses towards zero (set \code{keep.igamma = TRUE} for the
#'     diffuse inverse-gamma(0.001, 0.001) instead).
#'   \item JAGS used a uniform prior on q; JABBA3 uses a flat prior on log q (set \code{q.uniform = TRUE} to reproduce
#'     JABBA exactly).
#'   \item Plim (hockey-stick production) is not supported and is dropped.
#'   \item Auxiliary data are carried over with their type, SEs and the JABBA observation error settings
#'     (\code{auxiliary.sigma} -> estimated additional SD, \code{auxiliary.obsE} -> \code{aux.se.min}); JABBA's
#'     \code{auxiliary.lag} (mean exploitation over the previous years) has no JABBA3 equivalent and is reported
#'     (JABBA3 relates each observation to the same year). JABBA's uniform prior on qA becomes a flat prior on log qA.
#'   \item Catch-only models (\code{cpue = NULL}) become \code{jb3_data(catch, NULL)}; the status prior
#'     (\code{b.prior}) carries over.
#' }
#' @param jbinput output of \code{JABBA::build_jabba()}
#' @param q.uniform reproduce JABBA's uniform prior on q
#' @param keep.igamma use the diffuse inverse-gamma prior on the process variance instead of the FishLife prior
#' @param verbose print the conversion messages
#' @return list(data, priors, model), e.g. \code{x <- jbinput2jabba3(jb); fit <- jabba3(x$data, x$priors, x$model)}
#' @export
jbinput2jabba3 <- function(jbinput, q.uniform = FALSE, keep.igamma = FALSE, verbose = TRUE){
  s <- jbinput$settings; jd <- jbinput$jagsdata; dd <- jbinput$data
  msg <- function(...) if(verbose) message("jbinput2jabba3: ", ...)
  if(!is.null(jd$Plim) && jd$Plim > 0) msg("Plim (hockey-stick production) is not supported and was ignored")
  # data (catch-only and auxiliary data)
  co <- isTRUE(s$CatchOnly)
  aux <- if(isTRUE(s$Auxiliary)) dd$auxiliary else NULL
  atype <- s$auxiliary.type %||% "z"
  if(!is.null(aux)){
    msg("auxiliary data of type '", atype, "' carried over")
    if(!is.null(s$auxiliary.lag) && s$auxiliary.lag > 0 && atype %in% c("effort", "z", "ffmsy", "f"))
      msg("auxiliary.lag = ", s$auxiliary.lag, " has no JABBA3 equivalent: observations relate to the same year (lag 0)")
  }
  if(co) msg("catch-only model: no abundance index")
  data <- jb3_data(dd$catch, if(co) NULL else dd$cpue, if(isTRUE(s$SE.I) && !co) dd$se else NULL, stock = s$assessment %||% "stock",
                   auxiliary = aux, auxiliary.se = if(!is.null(aux)) dd$auxiliary.se else NULL, auxiliary.type = atype)
  # shape
  mt <- s$model.type %||% "Schaefer"
  shape <- switch(mt, Schaefer = "schaefer", Fox = "fox", "pella")
  bmsyk <- s$BmsyK %||% 0.4
  if(shape == "pella" && !is.null(jd$mu.m)) bmsyk <- jd$mu.m^(-1/(jd$mu.m - 1))
  # priors (JABBA jagsdata hold c(median, log.sd) for r, K and the lognormal psi)
  bk <- if(identical(s$psi.dist, "beta")){ a <- jd$psi.pr[1]; b <- jd$psi.pr[2]; mu <- a/(a + b); v <- a*b/((a + b)^2*(a + b + 1))
    msg("beta prior on psi converted to a lognormal prior on bkinit with the same mean and CV"); c(mu, cv2sd(sqrt(v)/mu))
  } else as.numeric(jd$psi.pr[1:2])
  m.pr <- if(identical(s$model.id, 4) || identical(mt, "Pella_m")) c(jd$mu.m, cv2sd(jd$m.CV %||% 0.3), 1) else c(NA, 0.3, 0)
  sdb <- if(!is.null(s$sigma.proc.prior)) as.numeric(s$sigma.proc.prior[1:2]) else if(keep.igamma) "igamma" else {
    msg("JABBA inverse-gamma prior on the process variance replaced by the JABBA3 FishLife prior on sdb"); "fishlife" }
  q <- NULL
  if(identical(s$q1.dist, "lnorm") && !is.null(jd$q1.pr)){
    nq <- jd$nq %||% 1; q <- cbind(c(jd$q1.pr[1], rep(1, nq - 1)), c(jd$q1.pr[2], rep(1, nq - 1)), c(1, rep(0, nq - 1)))
  }
  status <- NULL
  if(!is.null(jd$b.pr) && jd$b.pr[3] != 0){
    yrs <- data$years[which(jd$b.yr == 1)]
    type <- if(jd$b.pr[4] < 1) "BK" else if(jd$b.pr[4] > 1) "FFmsy" else "BBmsy"
    status <- data.frame(year = yrs, type = type, median = jd$b.pr[1], log.sd = jd$b.pr[2])
  }
  kd <- NULL
  if(!is.null(s$proc.k.prior)) kd <- c(s$proc.k.prior[1], s$proc.k.prior[2], 1)
  else if(!is.null(s$proc.floor)) kd <- c(s$proc.floor[2], 0.15, 0)
  priors <- jb3_priors(data, r = jd$r.pr[1:2], K = c(jd$K.pr[1], jd$K.pr[2], 1), bkinit = bk, sdb = sdb, q = q, m = m.pr,
                       kdyn = kd, status = status)
  # model
  fixed.sdb <- if(!is.null(s$proc.floor)) s$proc.floor[1] else if(is.numeric(s$sigma.proc)) s$sigma.proc else NULL
  if(!is.null(s$proc.floor)) msg("proc.floor (Tier 2) mapped to the dynamic sdb form with the floor level fixed")
  if(isTRUE(s$proc.mode == "dyn")) msg("proc.dyn (empirical ICES curve) has no JABBA3 equivalent; constant sdb used")
  se.min <- if(co) NULL else if(isTRUE(s$SE.I)) 0.001 else sqrt(min(jd$SE2, na.rm = TRUE))
  aux.se.min <- NULL; aux.sd <- NULL; aux.q <- NULL
  if(!is.null(aux)){
    if(is.null(dd$auxiliary.se) && !is.null(jd$A.SE2)) aux.se.min <- sqrt(min(jd$A.SE2, na.rm = TRUE))
    if(isFALSE(s$auxiliary.sigma)) aux.sd <- 1e-4
    aux.q <- rep(NA, ncol(aux) - 1)             # JABBA estimates qA for all auxiliary types
  }
  ebsb <- NULL; ebi <- NULL
  if(isTRUE(s$index.ebsb.on)){
    cu <- s$ebsb.setup$curves; ebsb <- lapply(seq_len(nrow(cu)), function(i) as.numeric(cu[i, 1:5])); names(ebsb) <- rownames(cu) %||% paste0("sel", seq_len(nrow(cu)))
    ebi <- as.integer(s$ebsb.setup$index)
    if(length(ebsb) > 1) msg("several EB/SB curves: removals use the first curve (multi-fleet catch shares not yet supported)")
  }
  model <- jb3_model(shape = shape, BmsyK = bmsyk, lag = s$proc.lag %||% 0,
                     sdb.form = if(!is.null(s$proc.floor)) "dyn" else "const", sdb.fixed = fixed.sdb,
                     xlim = s$proc.xlim %||% c(0.2, 3), biascor = isTRUE(s$proc.biascor),
                     start = if(isFALSE(s$proc.dev.all)) "index" else "all", se.min = se.min,
                     sets.q = jd$sets.q, sets.var = jd$sets.var, catch.error = isTRUE(s$add.catch.CV), sdc = s$catch.cv %||% 0.1,
                     ebsb = ebsb, ebsb.index = ebi, ebsb.removal = isTRUE(s$ebsb.removal) || is.null(ebsb),
                     q.uniform = q.uniform, P.bound = jd$P_bound %||% c(0.02, 1.3),
                     sdi.bound = jd$sigmaobs_bound %||% 1, sdi.est = !isFALSE(s$sigma.est), aux.q = aux.q, aux.sd = aux.sd, aux.se.min = aux.se.min)
  if(!q.uniform) msg("q has a flat prior on log q (JABBA/JAGS: uniform on q; use q.uniform = TRUE to reproduce)")
  list(data = data, priors = priors, model = model)
}

#' Convert a SPiCT input list to JABBA3
#'
#' Translates a SPiCT input list (\code{inp}, as passed to \code{spict::fit.spict()}) into JABBA3 \code{data},
#' \code{priors} and \code{model}. SPiCT works in continuous time; JABBA3 is annual, so catches are summed and index
#' observations within a year averaged (geometric mean) by calendar year (\code{floor(time)}); within-year index timing
#' is not used.
#'
#' Years: the model years span the catch and the indices, as in SPiCT. Index years before or after the catch series are
#' kept, and the catch in those years is estimated with an F random walk (\code{jb3_data(catch.gaps = "estimate")}).
#'
#' Observation error: the SPiCT form is used (\code{obs.error = "spict"}): the SD of an index observation is
#' \code{sdi[index] * stdevfacI}, with \code{stdevfacI} passed as the index se (multipliers of sdi, as in SPiCT). For
#' several observations in a year the multiplier of their geometric mean is \code{sqrt(sum(stdevfacI^2))/n}. Index
#' groups sharing sdi (\code{mapsdi}) and q (\code{mapq}) are kept (\code{sets.var}, \code{sets.q}). Catch error is on,
#' with the catch SD \code{sdc * stdevfacC} by year (\code{catch.se}), \code{sdc} from the \code{logsdc} prior (or
#' initial value).
#'
#' Priors: \code{logr}, \code{logK}, \code{logbkfrac} (-> bkinit, first model year), \code{logsdb} (-> sdb), \code{logn}
#' (-> shape m), \code{logsdi} (-> sdi; one prior per index when given as a list, as in SPiCT) and \code{logsdf}
#' (-> sdf, F random walk in years without catch) are carried over with their on/off switches (prior vectors
#' \code{c(mean, sd[, on])} on the log scale). A fixed n (\code{phases$logn = -1}) fixes the shape. SPiCT's r and
#' JABBA's r coincide for the Pella-Tomlinson form (MSY = rK / n^(n/(n-1))). The ratio priors \code{logalpha},
#' \code{logbeta} and the F process in years with catch have no JABBA3 equivalent and are reported.
#' @param inp SPiCT input list
#' @param stock stock name
#' @param verbose print conversion messages
#' @param spict.like TRUE (default): SPiCT's process structure as far as JABBA3 has it: random walk on F in all years
#'   (\code{frw.all = TRUE}; with catch known, the smoothness of F is what identifies the scale K in SPiCT), catch SD
#'   estimated when SPiCT has a \code{logsdc} prior (\code{sdc.est = TRUE}; without a prior the estimate tends to collapse
#'   to 0 in JABBA3, so sdc stays fixed), sdf without prior unless SPiCT has one, and the
#'   best of prior-median and SPiCT-like starting values (\code{inits = "best"}). FALSE: catch-conditioned JABBA model
#'   with the catch SD fixed at sdc (earlier behaviour).
#' @return list(data, priors, model)
#' @export
spictinp2jabba3 <- function(inp, stock = "stock", verbose = TRUE, spict.like = TRUE){
  msg <- function(...) if(verbose) message("spictinp2jabba3: ", ...)
  pr <- inp$priors %||% list(); ini <- inp$ini %||% list(); ph <- inp$phases %||% list()
  if(is.null(inp$obsC) || is.null(inp$timeC)) stop("inp needs obsC and timeC", call. = FALSE)
  # catch: annual sums; catch SD multipliers (stdevfacC) averaged by year
  cy <- floor(inp$timeC + 1e-8)
  cyrs <- sort(unique(cy))
  catch <- data.frame(year = cyrs, catch = as.numeric(tapply(inp$obsC, cy, sum)))
  if(any(inp$timeC %% 1 > 1e-8)) msg("sub-annual catches summed by calendar year")
  sfC <- if(is.null(inp$stdevfacC)) rep(1, length(cyrs)) else as.numeric(tapply(inp$stdevfacC, cy, mean))
  # indices: lists or a single vector
  obsI <- if(is.list(inp$obsI)) inp$obsI else list(inp$obsI)
  timeI <- if(is.list(inp$timeI)) inp$timeI else list(inp$timeI)
  sfI <- inp$stdevfacI; if(!is.null(sfI) && !is.list(sfI)) sfI <- list(sfI)
  if(is.null(sfI)) sfI <- lapply(obsI, function(x) rep(1, length(x)))
  nm <- names(obsI); if(is.null(nm) || any(nm == "")) nm <- paste0("I", seq_along(obsI))
  iy <- lapply(timeI, function(t) floor(t + 1e-8))
  yrs <- seq(min(c(cyrs, unlist(iy))), max(c(cyrs, unlist(iy))))
  I <- data.frame(year = yrs); S <- data.frame(year = yrs); seas <- FALSE
  for(i in seq_along(obsI)){
    ty <- iy[[i]]; ok <- is.finite(obsI[[i]]) & obsI[[i]] > 0
    if(any(duplicated(ty[ok]))) seas <- TRUE
    gm <- tapply(log(obsI[[i]][ok]), ty[ok], mean)
    sm <- tapply(sfI[[i]][ok], ty[ok], function(f) sqrt(sum(f^2))/length(f))
    k <- match(yrs, as.numeric(names(gm)))
    I[[nm[i]]] <- exp(as.numeric(gm))[k]; S[[nm[i]]] <- as.numeric(sm)[k]
  }
  if(seas) msg("several index observations per year averaged (geometric mean; stdevfacI combined as sqrt(sum(f^2))/n)")
  if(any(unlist(timeI) %% 1 > 1e-8)) msg("within-year index timing not used (annual model: indices refer to the start of the year)")
  gaps <- length(setdiff(yrs, cyrs)) > 0
  if(gaps) msg("model years ", min(yrs), "-", max(yrs), " span catch and indices; catch in ", length(setdiff(yrs, cyrs)),
               " years without catch (", .yr_ranges(setdiff(yrs, cyrs)), ") estimated with an F random walk")
  # priors: SPiCT c(mean, sd, on) on the log scale -> JABBA3 c(median, log.sd, on); lists (by index) -> one row per index
  lp1 <- function(x){ x <- as.numeric(x); if(length(x) == 2) x <- c(x, 1); c(exp(x[1]), x[2], x[3]) }
  lp <- function(p, default){ x <- pr[[p]]; if(is.null(x)) return(default)
    if(is.list(x)){ m <- do.call(rbind, lapply(x, lp1)); rownames(m) <- names(x); return(m) }
    lp1(x) }
  on <- function(x) !is.null(x) && any(if(is.matrix(x)) x[, 3] == 1 else x[3] == 1)
  r <- lp("logr", c(0.2, 0.5, 1)); K <- lp("logK", c(NA, cv2sd(1), 0)); bk <- lp("logbkfrac", c(0.9, 0.25, 0))
  sdb <- lp("logsdb", NULL); n <- lp("logn", NULL); sdf <- lp("logsdf", NULL)
  if(K[3] == 0) K <- c(NA, cv2sd(1), 0)            # SPiCT placeholder for an unused K prior: JABBA3 default start (8 x max catch)
  if(r[3] == 0) msg("no r prior in SPiCT: JABBA3 runs with r flat on the log scale (SPiCT benchmarks normally use one)")
  if(any(vapply(c("logalpha", "logbeta"), function(p) !is.null(pr[[p]]) && length(pr[[p]]) > 2 && pr[[p]][3] == 1, TRUE)))
    msg("alpha/beta ratio priors have no JABBA3 equivalent and were dropped")
  if(!on(sdb)){ sdb <- "fishlife"; msg("no sdb prior: JABBA3 FishLife default prior on sdb used") }
  # observation-variance and catchability groups
  ni <- length(obsI)
  sets.var <- if(!is.null(inp$mapsdi) && length(inp$mapsdi) == ni) as.integer(inp$mapsdi) else seq_len(ni)
  if(isTRUE(inp$onesdi)) sets.var <- rep(1L, ni)
  sets.q <- if(!is.null(inp$mapq) && length(inp$mapq) == ni) as.integer(inp$mapq) else seq_len(ni)
  sdi <- lp("logsdi", c(0.2, 0.5, 0))
  if(is.matrix(sdi)){
    if(nrow(sdi) != ni) stop("logsdi has ", nrow(sdi), " priors for ", ni, " indices", call. = FALSE)
    sdi <- sdi[!duplicated(sets.var), , drop = FALSE]          # one prior per sdi group (first index of the group)
  }
  n.fix <- !is.null(ph$logn) && ph$logn < 0
  if(n.fix){ m <- exp(as.numeric(ini$logn %||% log(2))); m.pr <- c(NA, 0.3, 0) }
  else if(!is.null(n) && n[3] == 1){ m <- n[1]; m.pr <- n }
  else if(on(pr$logngamma)){                     # SPiCT gamma prior on log n (shape, rate): lognormal on n with the same mean/sd of log n
    g <- as.numeric(pr$logngamma); mu <- g[1]/g[2]; sdl <- sqrt(g[1])/g[2]
    m <- exp(mu); m.pr <- c(exp(mu), sdl, 1)
    msg("logngamma prior (", signif(g[1], 3), ", ", signif(g[2], 3), ") carried as m ~ LN(", signif(m, 3), ", ", signif(sdl, 2), ")") }
  else { m <- 2; m.pr <- c(NA, 0.3, 0); msg("n neither fixed nor with a prior: Schaefer shape (n = 2) used") }
  shape <- if(abs(m - 2) < 1e-6) "schaefer" else "pella"
  bmsyk <- if(shape == "schaefer") 0.5 else m^(-1/(m - 1))
  sdc <- if(on(pr$logsdc)) exp(pr$logsdc[1]) else exp(as.numeric(ini$logsdc %||% log(0.1)))
  cse <- data.frame(year = cyrs, se = sdc*sfC)
  data <- jb3_data(catch, I, S, stock = stock, catch.se = cse, catch.gaps = if(gaps) "estimate" else "stop")
  sdc.fix <- !is.null(ph$logsdc) && ph$logsdc < 0
  sdc.pr <- if(on(pr$logsdc)) lp1(pr$logsdc) else c(0.1, 1, 0)
  priors <- jb3_priors(data, r = r, K = K, bkinit = bk, sdb = sdb, m = m.pr, sdi = sdi,
                       sdf = if(on(sdf)) sdf else if(spict.like) c(0.3, 0.5, 0) else c(0.3, 0.5, 1), sdc = sdc.pr)
  model <- jb3_model(shape = shape, BmsyK = bmsyk, obs.error = "spict", catch.error = TRUE, sdc = sdc,
                     sets.var = sets.var, sets.q = sets.q, frw.all = spict.like, sdc.est = spict.like && !sdc.fix && on(pr$logsdc),
                     inits = if(spict.like) "best" else "prior")
  if(spict.like) msg("SPiCT-like: random walk on F in all years (frw.all = TRUE; sdf ", if(on(sdf)) "with the SPiCT prior" else "without prior",
                     "); catch SD ", if(on(pr$logsdc) && !sdc.fix) "estimated with the SPiCT logsdc prior (sdc.est = TRUE; stdevfacC as multipliers)" else
                       paste0("fixed at ", signif(sdc, 2), " (no SPiCT logsdc prior: estimated without a prior, sdc tends to collapse to 0)"),
                     "; best of prior-median and SPiCT-like starts (inits = \"best\")")
  else msg("catch error: sdc = ", signif(sdc, 2), if(any(sfC != 1)) paste0(" x stdevfacC (catch SD ", paste(signif(range(sdc*sfC), 2), collapse = "-"), ")"),
      "; the F process is only modelled in years without catch")
  list(data = data, priors = priors, model = model)
}
