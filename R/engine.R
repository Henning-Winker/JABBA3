# RTMB engine: setup, objective, starting values

# Build the data list, parameters, map and random effects
jb3_setup <- function(data, priors, model, peels = 0, index.peels = 0){
  N <- length(data$years); nI <- ncol(data$I)
  if(peels > 0){                                  # retrospective: drop the last `peels` years of all data
    keep <- seq_len(N - peels)
    data$years <- data$years[keep]; data$catch <- data$catch[keep]
    data$I <- data$I[keep, , drop = FALSE]; data$SE <- data$SE[keep, , drop = FALSE]; N <- length(keep)
  }
  I <- data$I
  if(index.peels > 0) I[(N - index.peels + 1):N, ] <- NA   # hindcast: drop index values only
  mucatch <- mean(data$catch)
  ok <- which(!is.na(I), arr.ind = TRUE); ok <- ok[order(ok[, 1], ok[, 2]), , drop = FALSE]
  if(!nrow(ok)) stop("no index observations left", call. = FALSE)
  sets.q <- model$sets.q %||% seq_len(nI); sets.var <- model$sets.var %||% seq_len(nI)
  nq <- length(unique(sets.q)); nvar <- length(unique(sets.var))
  # base observation variance (JABBA 2 conventions) and SPiCT-type relative weights
  hasse <- colSums(!is.na(data$SE)) > 0
  semin <- model$se.min %||% NULL
  SE2 <- W <- matrix(1, N, nI)
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
  fi <- min(ok[, 1])
  dat <- list(
    N = N, TC = pmax(data$catch/mucatch, 1e-8), logIobs = log(I[ok]), obs_t = ok[, 1], obs_i = ok[, 2],
    SE2 = SE2, W = W, nI = nI, nq = nq, nvar = nvar, sets.q = as.integer(factor(sets.q)), sets.var = as.integer(factor(sets.var)),
    lagidx = pmax(seq_len(N) - model$lag, 1), lag0 = model$lag, stI = if(model$start == "index") fi else 1L,
    r.pr = priors$r, K.pr = c(K.pr[1]/mucatch, K.pr[2:3]), bk.pr = priors$bkinit, sdb.pr = priors$sdb,
    sdb.ig = priors$sdb.type == "igamma", q.pr = q.pr, m.pr = priors$m, mu.m = model$m, kd.pr = priors$kdyn,
    sdi.pr = priors$sdi,
    est.m = priors$m[3] == 1, est.sdb = is.null(model$sdb.fixed), sdb.fix = model$sdb.fixed %||% 0.1,
    rho.fix = model$rho, est.rho = priors$rho[3] == 1, rho.pr = priors$rho,
    dyn = model$sdb.form == "dyn", est.k = model$sdb.form == "dyn" && priors$kdyn[3] == 1, xlim = model$xlim,
    biascor = isTRUE(model$biascor), spict = model$obs.error == "spict",
    catch.err = isTRUE(model$catch.error), sdc = model$sdc, q.uniform = isTRUE(model$q.uniform),
    eb.on = eb.on, ebm = ebm, ebci = as.integer(ei), ebcc = as.integer(ec), eb.rem = eb.on && isTRUE(model$ebsb.removal),
    nst = if(is.null(st)) 0L else nrow(st), st.t = if(is.null(st)) 1L else match(st$year, data$years),
    st.type = if(is.null(st)) "BK" else as.character(st$type), st.pr = if(is.null(st)) matrix(1, 1, 2) else cbind(st$median, st$log.sd),
    P.bound = model$P.bound, K.bound = model$K.bound/mucatch, sdi.bound = model$sdi.bound
  )
  if(dat$nst > 0 && any(is.na(dat$st.t))) stop("status prior years outside the data years", call. = FALSE)
  par <- jb3_inits(dat)
  map <- list()
  if(!dat$est.m) map$logm <- factor(NA)
  if(!dat$est.sdb) map$logsdb <- factor(NA)
  if(!dat$est.k) map$tk <- factor(NA)
  if(!dat$est.rho) map$trho <- factor(NA)
  if(!dat$catch.err) map$cdev <- factor(rep(NA, N))
  random <- c("logP", if(dat$catch.err) "cdev")
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
  sdb0 <- if(dat$sdb.ig) 0.1 else dat$sdb.pr[1]
  k0 <- min(max((dat$kd.pr[1] + 0.5)/2.5, 0.05), 0.95)
  lapply(list(logK = log(K), logr = log(r), logq = lq, lbk = log(bk), logm = log(m),
       logsdb = log(if(dat$est.sdb) sdb0 else dat$sdb.fix), logsdi = rep(log(0.2), dat$nvar), tk = stats::qlogis(k0), trho = stats::qlogis(min(max((if(dat$est.rho) dat$rho.pr[1] else dat$rho.fix)/0.99, 0.01), 0.99)),
       logP = log(P), cdev = rep(0, N)), unname)
}

# Negative log posterior (RTMB)
jb3_nll <- function(par, dat){
  RTMB::getAll(par, dat, warn = FALSE)
  dnorm <- RTMB::dnorm
  logIobs <- RTMB::OBS(logIobs)
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
  else nll <- nll - sum(0.001*(-2*logsdi) - 0.001*exp(-2*logsdi))         # JABBA 2: gamma(0.001, 0.001) on 1/sdi^2
  if(dyn){
    if(est.k){ pk <- 1/(1 + exp(-tk)); kdyn <- -0.5 + 2.5*pk
      nll <- nll - dnorm(kdyn, kd.pr[1], kd.pr[2], TRUE) - log(2.5*pk*(1 - pk)) } else kdyn <- kd.pr[1]
  } else kdyn <- 0
  #--- catch
  if(catch.err){ estC <- TC*exp(cdev); nll <- nll - sum(dnorm(cdev, 0, sdc, TRUE)) } else estC <- TC
  #--- reference points and states
  BmsyK <- m^(-1/(m - 1)); Hmsy <- r/(m - 1)*(1 - 1/m)
  P <- exp(logP); BBmsy <- P/BmsyK
  rhof <- function(p, c) ebm[c, 1] + (ebm[c, 2] - ebm[c, 1])*(1 - exp(-ebm[c, 3]*(clamp(p, ebm[c, 4], ebm[c, 5], 1e-4) - ebm[c, 4])))/
    (1 - exp(-ebm[c, 3]*(ebm[c, 5] - ebm[c, 4])))
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
  if(est.rho){ prho <- 1/(1 + exp(-trho)); rho <- 0.99*prho
    nll <- nll - dnorm(rho, rho.pr[1], rho.pr[2], TRUE) - log(0.99*prho*(1 - prho)) } else rho <- rho.fix
  if(est.rho || rho.fix > 0){
    # AR(1) process deviations with marginal SD sdb_t: d_t ~ N(rho d_{t-1}, v_t (1 - rho^2)), first transition N(0, v_2)
    d <- logP - mu
    nll <- nll - dnorm(logP[1], mu[1], sqrt(v[1]), TRUE)
    if(N > 1) nll <- nll - dnorm(d[2], 0, sqrt(v[2]), TRUE)
    if(N > 2) nll <- nll - sum(dnorm(d[3:N], rho*d[2:(N - 1)], sqrt(v[3:N]*(1 - rho^2) + 1e-8), TRUE))
  } else nll <- nll - sum(dnorm(logP, mu, sqrt(v), TRUE))
  procdev <- logP - Pmean
  #--- soft penalties (JABBA 2)
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
  #--- observation equation
  ebi <- 0*logP[obs_t] + 1
  if(eb.on) for(i in seq_len(nI)){ w <- which(obs_i == i); if(length(w)) ebi[w] <- rhof(P[obs_t[w]], ebci[i]) }
  Imean <- logq[sets.q[obs_i]] + logP[obs_t] + logK + log(ebi)
  sdI <- if(spict) sdi[sets.var[obs_i]]*W[cbind(obs_t, obs_i)] else sqrt(SE2[cbind(obs_t, obs_i)] + sdi[sets.var[obs_i]]^2)
  nll <- nll - sum(dnorm(logIobs, Imean, sdI, TRUE))
  #--- outputs
  logB <- logP + logK; logBBmsy <- log(BBmsy); logFFmsy <- log(HHmsy); logF <- log(H)
  logBmsy <- log(BmsyK) + logK; logFmsy <- log(Hmsy); logMSY <- log(MSY); logbk <- lbk
  RTMB::ADREPORT(logB); RTMB::ADREPORT(logBBmsy); RTMB::ADREPORT(logFFmsy); RTMB::ADREPORT(logF); RTMB::ADREPORT(Imean)
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
  nll
}
