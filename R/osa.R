#' One-step-ahead residuals
#'
#' OSA (prediction) residuals of the abundance indices, as in SPiCT: residual = qnorm(P(I_t <= obs | I_1..I_t-1)) with
#' the fixed effects at the mode and the states integrated sequentially by a bootstrap particle filter (exact up to
#' Monte Carlo error for the nonlinear production model, including lag, dynamic process error and Selex). Under a
#' correct model they are iid N(0, 1). Ordinary residuals from the smoothed states are shrunk and autocorrelated by
#' construction, so residual tests use the OSA residuals by default.
#' @param fit a \code{jabba3} fit
#' @param np number of particles
#' @param seed random seed
#' @return numeric vector, in the order of the observations in \code{fit$index}
#' @export
jb3_osa <- function(fit, np = 5000, seed = 1){
  set.seed(seed)
  d <- fit$setup$dat; rep <- fit$report; N <- d$N
  K <- rep$K; r <- rep$r; m <- rep$m; BmsyK <- rep$BmsyK; sdb <- rep$sdb; kd <- rep$kdyn
  C <- rep$estC; lq <- log(rep$q); sdi <- rep$sdi
  rem <- function(p) if(d$eb.rem) .rho_curve(p, d$ebm[d$ebcc, ]) else 1
  sdf <- function(x) if(d$dyn) sdb*pmin(pmax(x, d$xlim[1]), d$xlim[2])^(-kd) else rep(sdb, length(x))
  v1 <- if(d$stI > 1) 1e-4 else sdf(rep$bk/BmsyK)[1]^2
  LP <- matrix(NA, np, N); LP[, 1] <- stats::rnorm(np, log(rep$bk), sqrt(v1))
  w <- rep(1/np, np); u <- rep(NA, length(d$logIobs))
  for(t in seq_len(N)){
    P <- exp(LP[, t])
    for(k in which(d$obs_t == t)){
      i <- d$obs_i[k]; y <- d$logIobs[k]
      mu <- lq[d$sets.q[i]] + LP[, t] + log(K) + (if(d$eb.on) log(.rho_curve(P, d$ebm[d$ebci[i], ])) else 0)
      sd <- if(d$spict) sdi[d$sets.var[i]]*d$W[t, i] else sqrt(d$SE2[t, i] + sdi[d$sets.var[i]]^2)
      u[k] <- sum(w*stats::pnorm(y, mu, sd))
      w <- w*stats::dnorm(y, mu, sd); w <- w/sum(w)
    }
    if(1/sum(w^2) < np/2){ id <- sample.int(np, np, TRUE, prob = w); LP <- LP[id, , drop = FALSE]; w <- rep(1/np, np); P <- exp(LP[, t]) }
    if(t < N){
      Pl <- exp(LP[, d$lagidx[t]])
      pm <- log(pmax(P + r/(m - 1)*Pl*(1 - Pl^(m - 1)) - C[t]/(K*rem(P)), 0.001))
      v <- if(t + 1 < d$stI) rep(1e-4, np) else sdf(P/BmsyK)^2
      LP[, t + 1] <- stats::rnorm(np, pm - if(d$biascor) 0.5*v else 0, sqrt(v))
    }
  }
  stats::qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6))
}
