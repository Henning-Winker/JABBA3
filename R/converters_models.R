#' Convert a JABBA 2 input object to JABBA3
#'
#' Translates a \code{JABBA::build_jabba()} object into the three JABBA3 inputs: \code{data}, \code{priors} and
#' \code{model}. Priors are carried over on the JABBA 2 scale (r, K, psi -> bkinit, shape, q1, status prior), as are the
#' lag, catch error, catchability and variance groups, the process-error start, the BioProcDyn options (lognormal
#' sigma.proc prior, fixed sigma.proc, depletion-dependent floor, bias correction) and the Selex EB/SB curves.
#'
#' Differences that cannot be carried over one to one are reported as messages:
#' \itemize{
#'   \item JABBA 2's inverse-gamma prior on the process variance is replaced by the JABBA3 default (FishLife lognormal
#'     prior on sdb), because the inverse-gamma mode collapses towards zero (set \code{keep.igamma = TRUE} for the
#'     diffuse inverse-gamma(0.001, 0.001) instead).
#'   \item JAGS used a uniform prior on q; JABBA3 uses a flat prior on log q (set \code{q.uniform = TRUE} to reproduce
#'     JABBA 2 exactly).
#'   \item Auxiliary data, catch-only mode and Plim are not yet supported and are dropped.
#' }
#' @param jbinput output of \code{JABBA::build_jabba()}
#' @param q.uniform reproduce JABBA 2's uniform prior on q
#' @param keep.igamma use the diffuse inverse-gamma prior on the process variance instead of the FishLife prior
#' @param verbose print the conversion messages
#' @return list(data, priors, model), e.g. \code{x <- jbinput2jabba3(jb); fit <- jabba3(x$data, x$priors, x$model)}
#' @export
jbinput2jabba3 <- function(jbinput, q.uniform = FALSE, keep.igamma = FALSE, verbose = TRUE){
  s <- jbinput$settings; jd <- jbinput$jagsdata; dd <- jbinput$data
  msg <- function(...) if(verbose) message("jbinput2jabba3: ", ...)
  if(isTRUE(s$CatchOnly)) stop("catch-only JABBA models are not supported by JABBA3 yet", call. = FALSE)
  if(isTRUE(s$Auxiliary)) msg("auxiliary data are not supported yet and were dropped")
  if(!is.null(jd$Plim) && jd$Plim > 0) msg("Plim (hockey-stick production) is not supported and was ignored")
  # data
  data <- jb3_data(dd$catch, dd$cpue, if(isTRUE(s$SE.I)) dd$se else NULL, stock = s$assessment %||% "stock")
  # shape
  mt <- s$model.type %||% "Schaefer"
  shape <- switch(mt, Schaefer = "schaefer", Fox = "fox", "pella")
  bmsyk <- s$BmsyK %||% 0.4
  if(shape == "pella" && !is.null(jd$mu.m)) bmsyk <- jd$mu.m^(-1/(jd$mu.m - 1))
  # priors (JABBA 2 jagsdata hold c(median, log.sd) for r, K and the lognormal psi)
  bk <- if(identical(s$psi.dist, "beta")){ a <- jd$psi.pr[1]; b <- jd$psi.pr[2]; mu <- a/(a + b); v <- a*b/((a + b)^2*(a + b + 1))
    msg("beta prior on psi converted to a lognormal prior on bkinit with the same mean and CV"); c(mu, cv2sd(sqrt(v)/mu))
  } else as.numeric(jd$psi.pr[1:2])
  m.pr <- if(identical(s$model.id, 4) || identical(mt, "Pella_m")) c(jd$mu.m, cv2sd(jd$m.CV %||% 0.3), 1) else c(NA, 0.3, 0)
  sdb <- if(!is.null(s$sigma.proc.prior)) as.numeric(s$sigma.proc.prior[1:2]) else if(keep.igamma) "igamma" else {
    msg("JABBA 2 inverse-gamma prior on the process variance replaced by the JABBA3 FishLife prior on sdb"); "fishlife" }
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
  se.min <- if(isTRUE(s$SE.I)) 0.001 else sqrt(min(jd$SE2, na.rm = TRUE))
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
                     sdi.bound = jd$sigmaobs_bound %||% 1)
  if(!q.uniform) msg("q has a flat prior on log q (JABBA 2/JAGS: uniform on q; use q.uniform = TRUE to reproduce)")
  list(data = data, priors = priors, model = model)
}

#' Convert a SPiCT input list to JABBA3
#'
#' Translates a SPiCT input list (\code{inp}, as passed to \code{spict::fit.spict()}) into JABBA3 \code{data},
#' \code{priors} and \code{model}. SPiCT works in continuous time; JABBA3 is annual, so catches are summed and index
#' observations within a year averaged (geometric mean) by calendar year (\code{floor(time)}). The SPiCT observation-error
#' form is used (\code{obs.error = "spict"}), with \code{stdevfacI} as relative year weights.
#'
#' Priors: \code{logr}, \code{logK}, \code{logbkfrac} (-> bkinit), \code{logsdb} (-> sdb), \code{logn} (-> shape m),
#' \code{logsdi} (-> sdi) are carried over with their on/off switches (prior vectors \code{c(mean, sd[, on])} on the log
#' scale). A fixed n (\code{phases$logn = -1}) fixes the shape. SPiCT's r and JABBA's r coincide for the Pella-Tomlinson
#' form (MSY = rK / n^(n/(n-1))). The ratio priors \code{logalpha}, \code{logbeta} and the F process have no JABBA3
#' equivalent and are reported. Catch error is on with \code{sdc} from the \code{logsdc} prior or initial value.
#' @param inp SPiCT input list
#' @param stock stock name
#' @param verbose print conversion messages
#' @return list(data, priors, model)
#' @export
spictinp2jabba3 <- function(inp, stock = "stock", verbose = TRUE){
  msg <- function(...) if(verbose) message("spictinp2jabba3: ", ...)
  pr <- inp$priors %||% list(); ini <- inp$ini %||% list(); ph <- inp$phases %||% list()
  # catch: annual sums
  if(is.null(inp$obsC) || is.null(inp$timeC)) stop("inp needs obsC and timeC", call. = FALSE)
  cy <- floor(inp$timeC + 1e-8)
  catch <- data.frame(year = sort(unique(cy)), catch = as.numeric(tapply(inp$obsC, cy, sum)))
  if(any(inp$timeC %% 1 > 1e-8)) msg("sub-annual catches summed by calendar year")
  # indices: lists or a single vector
  obsI <- if(is.list(inp$obsI)) inp$obsI else list(inp$obsI)
  timeI <- if(is.list(inp$timeI)) inp$timeI else list(inp$timeI)
  sfI <- inp$stdevfacI; if(!is.null(sfI) && !is.list(sfI)) sfI <- list(sfI)
  nm <- names(obsI); if(is.null(nm) || any(nm == "")) nm <- paste0("I", seq_along(obsI))
  yrs <- seq(min(catch$year), max(catch$year))
  I <- data.frame(year = yrs); S <- data.frame(year = yrs); seas <- FALSE
  for(i in seq_along(obsI)){
    ty <- floor(timeI[[i]] + 1e-8); if(any(duplicated(ty))) seas <- TRUE
    gm <- tapply(log(obsI[[i]]), ty, mean)
    I[[nm[i]]] <- exp(as.numeric(gm))[match(yrs, as.numeric(names(gm)))]
    S[[nm[i]]] <- if(is.null(sfI)) NA else as.numeric(tapply(sfI[[i]], ty, mean))[match(yrs, as.numeric(names(gm)))]
  }
  if(seas) msg("several index observations per year averaged (geometric mean)")
  outside <- setdiff(unlist(lapply(timeI, floor)), yrs)
  if(length(outside)) msg("index years outside the catch years dropped")
  data <- jb3_data(catch, I, if(is.null(sfI)) NULL else S, stock = stock)
  # priors: SPiCT c(mean, sd, on) on the log scale -> JABBA3 c(median, log.sd, on)
  lp <- function(nm, default){ x <- pr[[nm]]; if(is.null(x)) return(default); x <- as.numeric(x); if(length(x) == 2) x <- c(x, 1)
    c(exp(x[1]), x[2], x[3]) }
  r <- lp("logr", c(0.2, 0.5, 1)); K <- lp("logK", c(NA, cv2sd(1), 0)); bk <- lp("logbkfrac", c(0.9, 0.25, 0))
  sdb <- lp("logsdb", NULL); sdi <- lp("logsdi", c(0.2, 0.5, 0)); n <- lp("logn", NULL)
  if(r[3] == 0) msg("no r prior in SPiCT: JABBA3 runs with r flat on the log scale (SPiCT benchmarks normally use one)")
  if(any(vapply(c("logalpha", "logbeta"), function(p) !is.null(pr[[p]]) && length(pr[[p]]) > 2 && pr[[p]][3] == 1, TRUE)))
    msg("alpha/beta ratio priors have no JABBA3 equivalent and were dropped")
  if(is.null(sdb) || sdb[3] == 0){ sdb <- "fishlife"; msg("no sdb prior: JABBA3 FishLife default prior on sdb used") }
  n.fix <- !is.null(ph$logn) && ph$logn < 0
  if(n.fix){ m <- exp(as.numeric(ini$logn %||% log(2))); m.pr <- c(NA, 0.3, 0) }
  else if(!is.null(n) && n[3] == 1){ m <- n[1]; m.pr <- n }
  else { m <- 2; m.pr <- c(NA, 0.3, 0); msg("n neither fixed nor with a prior: Schaefer shape (n = 2) used") }
  shape <- if(abs(m - 2) < 1e-6) "schaefer" else "pella"
  bmsyk <- if(shape == "schaefer") 0.5 else m^(-1/(m - 1))
  priors <- jb3_priors(data, r = r, K = K, bkinit = bk, sdb = sdb, m = m.pr, sdi = sdi)
  sdc <- if(!is.null(pr$logsdc) && length(pr$logsdc) > 2 && pr$logsdc[3] == 1) exp(pr$logsdc[1]) else exp(as.numeric(ini$logsdc %||% log(0.1)))
  model <- jb3_model(shape = shape, BmsyK = bmsyk, obs.error = "spict", catch.error = TRUE, sdc = sdc)
  msg("the F process (sdf) is not modelled; catch enters with lognormal error (sdc = ", signif(sdc, 2), ")")
  list(data = data, priors = priors, model = model)
}
