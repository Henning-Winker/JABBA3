#' Residual diagnostics: JABBA set + SPiCT set
#'
#' One table per index, from the OSA residuals by default:
#' \itemize{
#'   \item JABBA set: n, RMSE of the log residuals, SDNR of the standardised residuals, runs test (p) with the share of
#'     residuals outside the three-sigma limits.
#'   \item SPiCT set (as \code{plotspict.diagnostic}): bias t-test, Ljung-Box test (lag 4) and Shapiro-Wilk normality test.
#' }
#' Under a correct model OSA residuals are iid N(0, 1); with \code{type = "ordinary"} the tests run on the classic
#' standardised residuals (smoothed states: shrunk and autocorrelated by construction, so the tests are biased).
#' Residual tests check the error assumptions, they have little power to detect structural misspecification (e.g. an
#' omitted production lag).
#' @param fit a \code{jabba3} fit
#' @param type "osa" (default) or "ordinary"
#' @return data.frame with one row per index and a \code{flag} column (any test with p < 0.05); the convergence verdict
#'   is in attribute "convergence"
#' @export
jb3_diags <- function(fit, type = c("osa", "ordinary")){
  type <- match.arg(type)
  ix <- fit$index
  if(type == "osa" && all(!is.finite(ix$osa))) stop("no OSA residuals: refit with osa = TRUE", call. = FALSE)
  res <- do.call(rbind, lapply(split(ix, factor(ix$index, unique(ix$index))), function(d){
    d <- d[order(d$year), ]
    z <- if(type == "osa") d$osa else d$std.residual
    z <- z[is.finite(z)]; n <- length(z)
    rt <- jb3_runs(if(type == "osa") z else d$residual)
    out3 <- if(type == "osa") mean(abs(z) > 3) else mean(d$residual < rt$sig3lim[1] | d$residual > rt$sig3lim[2])
    data.frame(index = d$index[1], n = n, RMSE = sqrt(mean(d$residual^2)), SDNR = stats::sd(z),
               runs.p = rt$p.runs, out.3sig = out3,
               bias.p = if(n > 2) stats::t.test(z)$p.value else NA,
               LjungBox.p = if(n > 5) stats::Box.test(z, lag = min(4, n - 2), type = "Ljung-Box")$p.value else NA,
               Shapiro.p = if(n > 3) stats::shapiro.test(z)$p.value else NA)
  }))
  res$flag <- apply(res[, c("runs.p", "bias.p", "LjungBox.p", "Shapiro.p")] < 0.05, 1, any, na.rm = TRUE)
  rownames(res) <- NULL
  attr(res, "type") <- type
  attr(res, "convergence") <- fit$convergence$tier
  res
}

#' Retrospective analysis
#'
#' Refits the model with the last 1..\code{n} years of all data removed and computes Mohn's rho for B, F, B/Bmsy and F/Fmsy.
#' @param fit a \code{jabba3} fit
#' @param n number of peels
#' @return list of class \code{jb3retro}: \code{fits} (reference first), \code{ts} (long), \code{rho}
#' @export
jb3_retro <- function(fit, n = 5){
  fits <- c(list(fit), lapply(seq_len(n), function(p) tryCatch(jabba3(fit$data, fit$priors, fit$model, engine = "rtmb",
                                                                          osa = FALSE, run = fit$run, peels = p), error = function(e) NULL)))
  ts <- do.call(rbind, lapply(seq_along(fits), function(k){ f <- fits[[k]]; if(is.null(f)) return(NULL)
    d <- f$ts[f$ts$quantity %in% c("B", "F", "B/Bmsy", "F/Fmsy"), ]; d$peel <- k - 1; d }))
  rho <- do.call(rbind, lapply(split(ts, ts$quantity), function(d){
    rel <- vapply(setdiff(unique(d$peel), 0), function(p){ ey <- max(d$year[d$peel == p])
      (d$est[d$peel == p & d$year == ey] - d$est[d$peel == 0 & d$year == ey])/d$est[d$peel == 0 & d$year == ey] }, 0)
    data.frame(quantity = d$quantity[1], rho = mean(rel)) }))
  structure(list(fits = fits, ts = ts, rho = rho, run = fit$run), class = "jb3retro")
}

#' Hindcast cross-validation
#'
#' Refits with the index observations of the last 1..\code{n} years removed (catch kept) and predicts each removed
#' observation one step ahead; MASE as in ss3diags (Winker and Kell), with the naive prediction = last observed value.
#' @param fit a \code{jabba3} fit
#' @param n number of peels
#' @param MAE.base.adj minimum naive error for MASE.adj
#' @return list of class \code{jb3hindcast}: \code{pred} (index, year, obs, hat, peel), \code{mase}
#' @export
jb3_hindcast <- function(fit, n = 5, MAE.base.adj = 0.1){
  ix0 <- fit$index; yrs <- fit$years; N <- length(yrs)
  pred <- do.call(rbind, lapply(seq_len(n), function(p){
    f <- tryCatch(jabba3(fit$data, fit$priors, fit$model, osa = FALSE, run = fit$run, index.peels = p), error = function(e) NULL)
    if(is.null(f)) return(NULL)
    d <- f$setup$dat; rp <- f$report
    # predicted index for all years and indices from the peeled fit
    out <- do.call(rbind, lapply(seq_len(d$nI), function(i){
      rho <- if(d$eb.on) .rho_curve(rp$P, d$ebm[d$ebci[i], ]) else 1
      data.frame(index = f$setup$index.names[i], year = yrs, hat = rp$q[d$sets.q[i]]*rp$P*rp$K*rho, peel = p)
    }))
    out[out$year <= yrs[N] - p + 1, ]
  }))
  obs <- ix0[, c("index", "year", "obs")]
  pred <- merge(pred, obs, by = c("index", "year"), all.x = TRUE)
  r0 <- data.frame(index = ix0$index, year = ix0$year, hat = ix0$hat, peel = 0, obs = ix0$obs)
  pred <- rbind(r0, pred[, names(r0)])
  ty <- yrs[N]; ev <- ty - c(0, seq_len(n))
  mase <- do.call(rbind, lapply(split(pred, pred$index), function(di){
    ref <- di[di$peel == 0, ]; ow <- ref[ref$year %in% ev, ]; ow <- ow[order(ow$year), ]
    pr <- do.call(rbind, lapply(seq_len(n), function(p){ yt <- ty - p + 1
      h <- di$hat[di$peel == p & di$year == yt]; o <- ref$obs[ref$year == yt]
      if(length(h) == 1 && length(o) == 1) data.frame(year = yt, pred.res = log(h) - log(o)) }))
    nv <- if(nrow(ow) > 1) data.frame(year = ow$year[-1], naive.res = diff(log(ow$obs)))
    r <- if(!is.null(pr) && !is.null(nv)) merge(pr, nv, by = "year") else NULL
    if(is.null(r) || !nrow(r)) return(data.frame(index = di$index[1], MASE = NA, MASE.adj = NA, MAE.PR = NA, MAE.base = NA, n.eval = 0))
    a <- mean(abs(r$pred.res)); b <- mean(abs(r$naive.res))
    data.frame(index = di$index[1], MASE = a/b, MASE.adj = a/max(b, MAE.base.adj), MAE.PR = a, MAE.base = b, n.eval = nrow(r))
  }))
  structure(list(pred = pred, mase = mase, run = fit$run, terminal = ty, n = n), class = "jb3hindcast")
}
