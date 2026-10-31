# RTMB engine: setup, objective, starting values

# Build the data list, parameters, map and random effects
jb3_setup <- function(data, priors, model, peels = 0, index.peels = 0, prior.only = FALSE, nc = FALSE){
  N <- length(data$years); nI <- ncol(data$I)
  if(peels > 0){                                  # retrospective: drop the last `peels` years of all data
    keep <- seq_len(N - peels)
    data$years <- data$years[keep]; data$catch <- data$catch[keep]
    data$I <- data$I[keep, , drop = FALSE]; data$SE <- data$SE[keep, , drop = FALSE]; N <- length(keep)
    if(!is.null(data$catch.se)) data$catch.se <- data$catch.se[keep]
    if(!is.null(data$A)){ data$A <- data$A[keep, , drop = FALSE]; data$A.SE <- data$A.SE[keep, , drop = FALSE] }
  }
  I <- data$I
  if(index.peels > 0) I[(N - index.peels + 1):N, ] <- NA   # hindcast: drop index values only
  mucatch <- mean(data$catch)
  ok <- which(!is.na(I), arr.ind = TRUE); ok <- ok[order(ok[, 1], ok[, 2]), , drop = FALSE]
  if(nI > 0 && !nrow(ok)) stop("no index observations left", call. = FALSE)
  sets.q <- model$sets.q %||% seq_len(nI); sets.var <- model$sets.var %||% seq_len(nI)
  if(nI == 0){ sets.q <- sets.var <- 1L }                 # catch-only: one dummy q and sdi (mapped off)
  nq <- length(unique(sets.q)); nvar <- length(unique(sets.var))
  # base observation variance (JABBA conventions) and SPiCT-type relative weights
  hasse <- colSums(!is.na(data$SE)) > 0
  semin <- model$se.min %||% NULL
  SE2 <- W <- matrix(1, N, max(nI, 1))
  for(i in seq_len(nI)){
    if(hasse[i]){ s <- data$SE[, i]; s[is.na(s)] <- model$se.na
      SE2[, i] <- s^2 + (semin %||% 0.001)^2
      W[, i] <- s/mean(s[!is.na(I[, i])])
    } else SE2[, i] <- (semin %||% 0.1)^2
  }
  # priors (K and q on the catch-scaled internal units)
  K.pr <- priors$K; if(is.na(K.pr[1])) K.pr[1] <- 8*max(data$catch)
  q.pr <- if(is.null(priors$q)) matrix(c(1, 1, 0), nq, 3, byrow = TRUE) else {
    x <- priors$q; if(nrow(x) == 1 && nq > 1) x <- x[rep(1, nq), , drop = FALSE]
    if(nrow(x) != nq) stop("q priors: one row per catchability group (", nq, ")", call. = FALSE); x[, 1] <- x[, 1]*mucatch; x }
  st <- priors$status
  # Selex curves
  eb.on <- !is.null(model$ebsb)
  if(eb.on){
    nc <- length(model$ebsb); ebm <- rbind(do.call(rbind, model$ebsb), SB = c(1, 1, 1, 0, 1))
    ei <- model$ebsb.index %||% rep(1, nI)
    if(is.character(ei)) ei <- match(ei, names(model$ebsb), nomatch = 0)
    ei <- rep(ei, length.out = nI); ei[ei == 0] <- nc + 1
    ec <- model$ebsb.catch; if(is.character(ec)) ec <- match(ec, names(model$ebsb))
  } else { ebm <- matrix(c(1, 1, 1, 0, 1), 1); ei <- rep(1, nI); ec <- 1 }
  fi <- if(nrow(ok)) min(ok[, 1]) else 1L
  ax <- .aux_setup(data, model, priors, N)
  dat <- list(
    N = N, TC = pmax(data$catch/mucatch, 1e-8), logIobs = log(I[ok]), obs_t = ok[, 1], obs_i = ok[, 2],
    SE2 = SE2, W = W, nI = nI, nq = nq, nvar = nvar, sets.q = as.integer(factor(sets.q)), sets.var = as.integer(factor(sets.var)),
    lagidx = pmax(seq_len(N) - model$lag, 1), lag0 = model$lag, stI = if(model$start == "index") fi else 1L,
    r.pr = priors$r, K.pr = c(K.pr[1]/mucatch, K.pr[2:3]), bk.pr = priors$bkinit, sdb.pr = priors$sdb,
    sdb.ig = priors$sdb.type == "igamma", q.pr = q.pr, m.pr = priors$m, mu.m = model$m, kd.pr = priors$kdyn,
    sdi.pr = priors$sdi,
    est.m = priors$m[3] == 1, est.sdb = is.null(model$sdb.fixed), sdb.fix = model$sdb.fixed %||% 0.1,
    rw = isTRUE(model$rw), rw.scale = model$rw.scale %||% 1,
    rho.fix = model$rho, est.rho = priors$rho[3] == 1, rho.pr = priors$rho,
    dyn = model$sdb.form == "dyn", est.k = model$sdb.form == "dyn" && priors$kdyn[3] == 1, xlim = model$xlim,
    biascor = isTRUE(model$biascor), spict = model$obs.error == "spict",
    catch.err = isTRUE(model$catch.error), sdc = .sdc_vec(model$sdc, data$catch.se, N), q.uniform = isTRUE(model$q.uniform),
    eb.on = eb.on, ebm = ebm, ebci = as.integer(ei), ebcc = as.integer(ec), eb.rem = eb.on && isTRUE(model$ebsb.removal),
    nst = if(is.null(st)) 0L else nrow(st), st.t = if(is.null(st)) 1L else match(st$year, data$years),
    st.type = if(is.null(st)) "BK" else as.character(st$type), st.pr = if(is.null(st)) matrix(1, 1, 2) else cbind(st$median, st$log.sd),
    P.bound = model$P.bound, K.bound = model$K.bound/mucatch, sdi.bound = model$sdi.bound,
    prior.only = isTRUE(prior.only), nc = isTRUE(nc)
  )
  dat <- c(dat, ax)
  if(dat$nst > 0){
    if(any(st$year[!is.na(st$year)] > max(data$years) + peels)) stop("status prior years outside the data years", call. = FALSE)
    dat$st.t[is.na(st$year)] <- N                          # NA year: terminal year of the fit
    if(any(is.na(dat$st.t))){
      if(peels > 0){ k <- !is.na(dat$st.t); dat$nst <- sum(k); dat$st.t <- if(any(k)) dat$st.t[k] else 1L
        dat$st.type <- if(any(k)) dat$st.type[k] else "BK"; dat$st.pr <- if(any(k)) dat$st.pr[k, , drop = FALSE] else matrix(1, 1, 2) }
      else stop("status prior years outside the data years", call. = FALSE)
    }
  }
  par <- jb3_inits(dat)
  map <- list()
  if(!dat$est.m) map$logm <- factor(NA)
  if(!dat$est.sdb) map$logsdb <- factor(NA)
  if(!dat$est.k) map$tk <- factor(NA)
  if(!dat$est.rho) map$trho <- factor(NA)
  if(!dat$catch.err) map$cdev <- factor(rep(NA, N))
  if(nI == 0 || dat$prior.only){ map$logq <- factor(rep(NA, nq)); map$logsdi <- factor(rep(NA, nvar)) }
  map$logqA <- factor(ifelse(dat$est.qA & !dat$prior.only, seq_along(dat$est.qA), NA))
  map$logsdA <- factor(ifelse(dat$est.sdA & !dat$prior.only, seq_along(dat$est.sdA), NA))
  if(!dat$est.creep) map$creep <- factor(NA)
  if(dat$nc){ par$zP <- rep(0, N); par$logP <- NULL }
  random <- c(if(dat$nc) "zP" else "logP", if(dat$catch.err) "cdev")
  list(dat = dat, par = par, map = map, random = random, mucatch = mucatch, years = data$years, obs = ok,
       index.names = colnames(data$I), data = data)
}

# starting values from a deterministic trajectory at the prior medians
jb3_inits <- function(dat, jitter = 0){
  N <- dat$N; m <- dat$mu.m
  K <- dat$K.pr[1]*exp(jitter*stats::rnorm(1, 0, 0.3)); r <- dat$r.pr[1]*exp(jitter*stats::rnorm(1, 0, 0.2))
  bk <- min(max(dat$bk.pr[1], 0.05), 1.2)
  rem <- function(p) if(dat$eb.rem) .rho_curve(p, dat$ebm[dat$ebcc, ]) else 1
  for(k in 1:15){
    P <- numeric(N); P[1] <- bk
    for(t in 2:N){ Pl <- P[dat$lagidx[t - 1]]; P[t] <- max(P[t-1] + r/(m - 1)*Pl*(1 - Pl^(m - 1)) - dat$TC[t-1]/(K*rem(P[t-1])), 0.02) }
    if(mean(P <= 0.021) < 0.1) break
    K <- K*1.5
  }
  rho <- function(i) if(dat$eb.on) .rho_curve(P[dat$obs_t], dat$ebm[dat$ebci[i], ]) else 1
  lq <- vapply(seq_len(dat$nq), function(j){
    w <- dat$sets.q[dat$obs_i] == j; if(!any(w)) return(log(1e-3))
    mean(dat$logIobs[w] - log(P[dat$obs_t[w]]*K*(if(dat$eb.on) .rho_curve(P[dat$obs_t[w]], dat$ebm[dat$ebci[dat$obs_i[w][1]], ]) else 1)))
  }, 0)
  H0 <- dat$TC/(K*P)
  lqA <- vapply(seq_len(max(dat$nA, 1)), function(i){
    if(dat$nA == 0 || !dat$est.qA[i]) return(log(dat$qA.fix[i] %||% 1))
    w <- dat$aobs_i == i; tt <- dat$aobs_t[w]
    x <- .aux_x(dat$atype[i], H0[tt], P[tt], P[tt]/dat$BmsyK0, H0[tt]/(r/(m - 1)*(1 - 1/m)))
    if(dat$atype[i] == 1) mean(log(x) - dat$logAobs[w] - (tt - dat$at0[i])*log(1 + dat$creep.fix)) else mean(dat$logAobs[w] - log(x))
  }, 0)
  sdb0 <- if(dat$sdb.ig) 0.1 else dat$sdb.pr[1]
  k0 <- min(max((dat$kd.pr[1] + 0.5)/2.5, 0.05), 0.95)
  lapply(list(logK = log(K), logr = log(r), logq = lq, lbk = log(bk), logm = log(m),
       logsdb = log(if(dat$est.sdb) sdb0 else dat$sdb.fix), logsdi = rep(log(0.2), dat$nvar), tk = stats::qlogis(k0), trho = stats::qlogis(min(max((if(dat$est.rho) dat$rho.pr[1] else dat$rho.fix)/0.99, 0.01), 0.99)),
       logP = log(P), cdev = rep(0, N), logqA = lqA, logsdA = log(dat$sdA.fix), creep = dat$creep.fix), unname)
}

# Negative log posterior (RTMB)
jb3_nll <- function(par, dat){
  RTMB::getAll(par, dat, warn = FALSE)
  dnorm <- RTMB::dnorm
  if(length(logIobs)) logIobs <- RTMB::OBS(logIobs)
  smax <- .smax; smin <- .smin; clamp <- .clamp; lclamp <- .lclamp
  nll <- 0
  #--- parameters and priors (log scale; lognormal priors = normal on the log)
  K <- exp(logK); r <- exp(logr); q <- exp(logq); bk <- exp(lbk)
  if(r.pr[3] == 1) nll <- nll - dnorm(logr, log(r.pr[1]), r.pr[2], TRUE)
  if(K.pr[3] == 1) nll <- nll - dnorm(logK, log(K.pr[1]), K.pr[2], TRUE)
  if(bk.pr[3] == 1) nll <- nll - dnorm(lbk, log(bk.pr[1]), bk.pr[2], TRUE)
  for(j in seq_len(nq)) if(q.pr[j, 3] == 1) nll <- nll - dnorm(logq[j], log(q.pr[j, 1]), q.pr[j, 2], TRUE)
  if(q.uniform) nll <- nll - sum(logq)
  if(est.m){ m <- exp(logm); nll <- nll - dnorm(logm, log(m.pr[1]), m.pr[2], TRUE) } else m <- mu.m
  sdb <- exp(logsdb)
  if(est.sdb){
    if(sdb.ig) nll <- nll - (0.001*(-2*logsdb) - 0.001*exp(-2*logsdb))       # gamma(0.001, 0.001) on 1/sdb^2 + Jacobian
    else if(sdb.pr[3] == 1) nll <- nll - dnorm(logsdb, log(sdb.pr[1]), sdb.pr[2], TRUE)
  }
  sdi <- exp(logsdi)
  if(spict){ if(sdi.pr[3] == 1) nll <- nll - sum(dnorm(logsdi, log(sdi.pr[1]), sdi.pr[2], TRUE)) }
  else nll <- nll - sum(0.001*(-2*logsdi) - 0.001*exp(-2*logsdi))         # JABBA: gamma(0.001, 0.001) on 1/sdi^2
  if(dyn){
    if(est.k){ pk <- 1/(1 + exp(-tk)); kdyn <- -0.5 + 2.5*pk
      nll <- nll - dnorm(kdyn, kd.pr[1], kd.pr[2], TRUE) - log(2.5*pk*(1 - pk)) } else kdyn <- kd.pr[1]
  } else kdyn <- 0
  #--- catch
  if(catch.err){ estC <- TC*exp(cdev); nll <- nll - sum(dnorm(cdev, 0, sdc, TRUE)) } else estC <- TC
  #--- reference points and states
  BmsyK <- m^(-1/(m - 1)); Hmsy <- r/(m - 1)*(1 - 1/m)
  rhof <- function(p, c) ebm[c, 1] + (ebm[c, 2] - ebm[c, 1])*(1 - exp(-ebm[c, 3]*(clamp(p, ebm[c, 4], ebm[c, 5], 1e-4) - ebm[c, 4])))/
    (1 - exp(-ebm[c, 3]*(ebm[c, 5] - ebm[c, 4])))
  if(est.rho){ prho <- 1/(1 + exp(-trho)); rho <- 0.99*prho
    nll <- nll - dnorm(rho, rho.pr[1], rho.pr[2], TRUE) - log(0.99*prho*(1 - prho)) } else rho <- rho.fix
  if(nc){
    # non-centred: standard normal innovations zP; the states follow from the process equation (same model as centred)
    sdf1 <- if(dyn) sdb*lclamp(bk/BmsyK, xlim[1], xlim[2])^(-kdyn) else sdb
    logP <- 0*zP
    logP[1] <- lbk + sqrt(if(stI > 1) 1e-4 else sdf1^2)*zP[1]
    dprev <- 0*lbk
    if(N > 1) for(t in 2:N){
      Pp <- exp(logP[t - 1]); Pl <- exp(logP[lagidx[t - 1]])
      remt <- if(eb.rem) rhof(Pp, ebcc) else 1
      pm <- log(smax(Pp + r/(m - 1)*Pl*(1 - Pl^(m - 1)) - estC[t - 1]/(K*remt), 0.001, 1e-7))
      sdp <- if(dyn) sdb*lclamp(Pp/BmsyK, xlim[1], xlim[2])^(-kdyn) else sdb
      vt <- if(t < stI) 1e-4 + 0*sdb else sdp^2
      mut <- if(biascor) pm - 0.5*vt else pm
      dt <- if(t == 2 || (!rw && !(est.rho || rho.fix > 0))) sqrt(vt)*zP[t] else
        if(rw) dprev + (rw.scale*sqrt(vt) + 1e-8)*zP[t] else rho*dprev + sqrt(vt*(1 - rho^2) + 1e-8)*zP[t]
      logP[t] <- mut + dt; dprev <- dt
    }
    nll <- nll - sum(dnorm(zP, 0, 1, TRUE))
  }
  P <- exp(logP); BBmsy <- P/BmsyK
  if(eb.rem){ rem <- rhof(P, ebcc); rem.msy <- rhof(BmsyK, ebcc) } else { rem <- 1 + 0*P; rem.msy <- 1 }
  MSY <- BmsyK*K*Hmsy*rem.msy
  #--- process variance (v[t] applies to the transition into year t)
  if(dyn){ sdt <- sdb*lclamp(BBmsy, xlim[1], xlim[2])^(-kdyn); sd1 <- sdb*lclamp(bk/BmsyK, xlim[1], xlim[2])^(-kdyn) }
  else { sdt <- sdb + 0*P; sd1 <- sdb }
  v <- 0*P
  v[1] <- if(stI > 1) 1e-4 else sd1^2
  if(N > 1) v[2:N] <- sdt[1:(N - 1)]^2
  tt <- seq_len(N); if(stI > 1) v[tt < stI] <- 1e-4
  #--- process equation
  Pl <- P[lagidx[-N]]
  Pmean <- c(lbk, log(smax(P[-N] + r/(m - 1)*Pl*(1 - Pl^(m - 1)) - estC[-N]/(K*rem[-N]), 0.001, 1e-7)))
  mu <- Pmean; if(biascor) mu[-1] <- Pmean[-1] - 0.5*v[-1]
  if(nc){
    # process likelihood already on zP
  } else if(rw){
    # random walk in the process deviations: d_t ~ N(d_{t-1}, (rw.scale sdb_t)^2), first transition N(0, v_2)
    d <- logP - mu
    nll <- nll - dnorm(logP[1], mu[1], sqrt(v[1]), TRUE)
    if(N > 1) nll <- nll - dnorm(d[2], 0, sqrt(v[2]), TRUE)
    if(N > 2) nll <- nll - sum(dnorm(d[3:N], d[2:(N - 1)], rw.scale*sqrt(v[3:N]) + 1e-8, TRUE))
  } else if(est.rho || rho.fix > 0){
    # AR(1) process deviations with marginal SD sdb_t: d_t ~ N(rho d_{t-1}, v_t (1 - rho^2)), first transition N(0, v_2)
    d <- logP - mu
    nll <- nll - dnorm(logP[1], mu[1], sqrt(v[1]), TRUE)
    if(N > 1) nll <- nll - dnorm(d[2], 0, sqrt(v[2]), TRUE)
    if(N > 2) nll <- nll - sum(dnorm(d[3:N], rho*d[2:(N - 1)], sqrt(v[3:N]*(1 - rho^2) + 1e-8), TRUE))
  } else nll <- nll - sum(dnorm(logP, mu, sqrt(v), TRUE))
  procdev <- logP - Pmean
  #--- soft penalties (JABBA)
  penB <- smin(logP - log(P.bound[1]), 0, 1e-4) + smax(logP - log(P.bound[2]), 0, 1e-4)
  nll <- nll + 500*sum(penB^2)
  penK <- smin(logK - log(K.bound[1]), 0, 1e-4) + smax(logK - log(K.bound[2]), 0, 1e-4)
  nll <- nll + 500*penK^2
  if(!spict) nll <- nll + 500*sum(smax(logsdi - log(sdi.bound), 0, 1e-4)^2)
  H <- estC/(K*P*rem); HHmsy <- H/Hmsy
  #--- status priors
  if(nst > 0) for(s in seq_len(nst)){
    x <- if(st.type[s] == "BK") P[st.t[s]] else if(st.type[s] == "BBmsy") BBmsy[st.t[s]] else HHmsy[st.t[s]]
    nll <- nll - dnorm(log(x), log(st.pr[s, 1]), st.pr[s, 2], TRUE)
  }
  #--- observation equation (indices)
  nobs <- length(obs_t)
  if(nobs > 0){
    ebi <- 0*logP[obs_t] + 1
    if(eb.on) for(i in seq_len(nI)){ w <- which(obs_i == i); if(length(w)) ebi[w] <- rhof(P[obs_t[w]], ebci[i]) }
    Imean <- logq[sets.q[obs_i]] + logP[obs_t] + logK + log(ebi)
    sdI <- if(spict) sdi[sets.var[obs_i]]*W[cbind(obs_t, obs_i)] else sqrt(SE2[cbind(obs_t, obs_i)] + sdi[sets.var[obs_i]]^2)
    if(!prior.only) nll <- nll - sum(dnorm(logIobs, Imean, sdI, TRUE))
  } else { Imean <- numeric(0); sdI <- numeric(0) }
  #--- auxiliary data (effort, F, F/Fmsy, B/Bmsy, B/K, Z)
  sdA <- exp(logsdA)
  if(nA > 0){
    Amean <- 0*logP[aobs_t]
    for(i in seq_len(nA)){
      w <- which(aobs_i == i); tt <- aobs_t[w]
      x <- if(atype[i] == 1) H[tt] else if(atype[i] == 2) -log(1 - smin(H[tt], 0.99, 1e-4)) else if(atype[i] == 3) HHmsy[tt] else
        if(atype[i] == 4) BBmsy[tt] else if(atype[i] == 5) P[tt] else smax(1 - (1 - Hmsy)*(1 - H[tt]), 0.001, 1e-6)
      Amean[w] <- if(atype[i] == 1) log(x) - logqA[i] - (tt - at0[i])*log(1 + creep) else logqA[i] + log(x)   # effort: E = F / (q (1 + creep)^(t - t0))
    }
    sdAo <- sqrt(ASE2 + sdA[aobs_i]^2)
    if(est.creep) nll <- nll - dnorm(creep, creep.pr[1], creep.pr[2], TRUE)
    if(!prior.only){
      nll <- nll - sum(dnorm(logAobs, Amean, sdAo, TRUE))
      for(i in seq_len(nA)) if(est.sdA[i]) nll <- nll - (0.001*(-2*logsdA[i]) - 0.001*exp(-2*logsdA[i])) +
        500*smax(logsdA[i] - log(sdi.bound), 0, 1e-4)^2
    }
  } else { Amean <- numeric(0); sdAo <- numeric(0) }
  #--- outputs
  logB <- logP + logK; logBBmsy <- log(BBmsy); logFFmsy <- log(HHmsy); logF <- log(H)
  logBmsy <- log(BmsyK) + logK; logFmsy <- log(Hmsy); logMSY <- log(MSY); logbk <- lbk
  RTMB::ADREPORT(logB); RTMB::ADREPORT(logBBmsy); RTMB::ADREPORT(logFFmsy); RTMB::ADREPORT(logF)
  if(nobs > 0) RTMB::ADREPORT(Imean)
  if(nA > 0){ RTMB::ADREPORT(Amean); RTMB::ADREPORT(logqA); RTMB::ADREPORT(logsdA) }
  if(est.creep) RTMB::ADREPORT(creep)
  if(nc) RTMB::ADREPORT(logP)
  RTMB::ADREPORT(logK); RTMB::ADREPORT(logr); RTMB::ADREPORT(logBmsy); RTMB::ADREPORT(logFmsy); RTMB::ADREPORT(logMSY)
  RTMB::ADREPORT(logbk); RTMB::ADREPORT(logsdb); RTMB::ADREPORT(logsdi); RTMB::ADREPORT(logq); RTMB::ADREPORT(procdev)
  if(est.m) RTMB::ADREPORT(logm)
  if(est.k) RTMB::ADREPORT(kdyn)
  if(est.rho) RTMB::ADREPORT(rho)
  if(catch.err) RTMB::ADREPORT(estC)
  RTMB::REPORT(K); RTMB::REPORT(r); RTMB::REPORT(q); RTMB::REPORT(bk); RTMB::REPORT(m); RTMB::REPORT(sdb); RTMB::REPORT(sdi)
  RTMB::REPORT(kdyn); RTMB::REPORT(rho); RTMB::REPORT(BmsyK); RTMB::REPORT(Hmsy); RTMB::REPORT(MSY); RTMB::REPORT(P); RTMB::REPORT(BBmsy)
  RTMB::REPORT(H); RTMB::REPORT(HHmsy); RTMB::REPORT(estC); RTMB::REPORT(Pmean); RTMB::REPORT(procdev); RTMB::REPORT(sdt)
  RTMB::REPORT(v); RTMB::REPORT(rem); RTMB::REPORT(rem.msy); RTMB::REPORT(Imean); RTMB::REPORT(sdI)
  qA <- exp(logqA)
  RTMB::REPORT(Amean); RTMB::REPORT(sdAo); RTMB::REPORT(qA); RTMB::REPORT(sdA); RTMB::REPORT(creep); RTMB::REPORT(mu); RTMB::REPORT(logP)
  nll
}

# year-specific catch SD (log scale): data$catch.se (from jb3_data(catch.se =), NA = model sdc) or a vector sdc in the model
.sdc_vec <- function(sdc, cse, N){
  v <- if(length(sdc) == N) sdc else rep(sdc[1], N)
  if(!is.null(cse)){ cse <- rep(cse, length.out = N); v[is.finite(cse) & cse > 0] <- cse[is.finite(cse) & cse > 0] }
  pmax(v, 1e-4)
}

# auxiliary data setup: observation vectors, types, catchability / SD / creep switches
.aux_setup <- function(data, model, priors, N){
  A <- data$A; tcode <- c(effort = 1L, f = 2L, ffmsy = 3L, bbmsy = 4L, bk = 5L, z = 6L)
  cr <- model$effort.creep %||% 0; est.creep <- identical(cr, "estimate")
  creep.pr <- priors$creep %||% c(0.02, 0.02)
  BmsyK0 <- model$BmsyK
  base <- list(nA = 0L, logAobs = numeric(0), aobs_t = integer(0), aobs_i = integer(0), ASE2 = numeric(0), atype = 1L,
               at0 = 1L, est.qA = FALSE, qA.fix = 1, est.sdA = FALSE, sdA.fix = 0.2, creep.fix = if(est.creep) as.numeric(creep.pr[1]) else as.numeric(cr),
               est.creep = FALSE, creep.pr = as.numeric(creep.pr[1:2]), BmsyK0 = BmsyK0,
               aux.names = character(0))
  if(is.null(A) || !ncol(A)) return(base)
  nA <- ncol(A); tp <- tcode[data$A.type]
  ok <- which(!is.na(A), arr.ind = TRUE); ok <- ok[order(ok[, 2], ok[, 1]), , drop = FALSE]
  hasse <- colSums(!is.na(data$A.SE)) > 0
  semin <- model$aux.se.min
  ase2 <- vapply(seq_len(nrow(ok)), function(k){ t <- ok[k, 1]; i <- ok[k, 2]
    if(hasse[i]){ s <- data$A.SE[t, i]; if(is.na(s)) s <- model$se.na; s^2 + (semin %||% 0.001)^2 } else (semin %||% 0.1)^2 }, 0)
  aq <- model$aux.q
  qfix <- if(is.null(aq)) ifelse(tp %in% c(1L, 6L), NA, 1) else rep(as.numeric(aq), length.out = nA)
  asd <- model$aux.sd
  sdfix <- if(is.null(asd)) rep(NA_real_, nA) else rep(as.numeric(asd), length.out = nA)
  base$nA <- nA; base$logAobs <- log(A[ok]); base$aobs_t <- as.integer(ok[, 1]); base$aobs_i <- as.integer(ok[, 2])
  base$ASE2 <- ase2; base$atype <- as.integer(tp)
  base$at0 <- vapply(seq_len(nA), function(i) as.integer(min(which(!is.na(A[, i])))), 1L)
  base$est.qA <- is.na(qfix); base$qA.fix <- ifelse(is.na(qfix), 1, qfix)
  base$est.sdA <- is.na(sdfix); base$sdA.fix <- ifelse(is.na(sdfix), 0.2, pmax(sdfix, 1e-4))
  base$est.creep <- est.creep && any(tp == 1L)
  base$aux.names <- colnames(A)
  base
}

.aux_x <- function(type, H, P, BBmsy, HHmsy) switch(type, H, -log(1 - pmin(H, 0.99)), HHmsy, BBmsy, P, pmax(1 - (1 - 0.1)*(1 - H), 0.001))
